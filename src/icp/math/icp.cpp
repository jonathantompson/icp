#include <random>
#include <stdexcept>
#include <iostream>
#include "icp/math/icp.h"
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include "icp/math/icp_eigen_data.h"
#include <nanoflann/nanoflann.hpp>
#include "icp/math/bfgs.h"
#include "icp/math/math_base.h"

#define SAFE_DELETE(x) if (x != NULL) { delete x; x = NULL; }
#define SAFE_DELETE_ARR(x) if (x != NULL) { delete[] x; x = NULL; }

using Eigen::MatrixXf;
using std::cout;
using std::endl;
using std::runtime_error;
using namespace nanoflann;

// nanoflann needs a DataSet Adaptor class to calculate distances
template <typename T>
struct PointCloud {
  PointCloud(T* points, const uint32_t npoints) : points(points), 
    npoints(npoints) { }
  T* points;
  uint32_t npoints;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return npoints; }

	// Returns the distance between the vector "p1[0:size-1]" and the data point
  // with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T* p1, const size_t idx_p2, size_t size) const {
		const T d0=p1[0] - points[idx_p2 * 3 + 0];
		const T d1=p1[1] - points[idx_p2 * 3 + 1];
		const T d2=p1[2] - points[idx_p2 * 3 + 2];
		return d0 * d0 + d1 * d1 + d2 * d2;  // L2 squared
	}

	// Returns the dim'th component of the idx'th point in the class
	inline T kdtree_get_pt(const size_t idx, int dim) const {
    return points[idx * 3 + dim];
	}

	// Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }
};

namespace icp {
namespace math {

  template <typename T> icp::data_str::Vector<double> ICP<T>::cur_Q_;
  template <typename T> icp::data_str::Vector<double> ICP<T>::cur_D_;
  template <typename T> icp::data_str::Vector<double> ICP<T>::cur_weights_;
  template <typename T> Double4x4 ICP<T>::cur_mat_;

  template <typename T>
  ICP<T>::ICP() {
    bfgs_ = NULL;
    edata_ = NULL;

    edata_ = new ICPEigenData();
    bfgs_ = new BFGS<double>(ICP_BFGS_NUM_COEFFS);
    verbose = ICP_DEFAULT_VERBOSE;
    num_iterations = ICP_DEFAULT_ITERATIONS;
    cos_normal_threshold = (T)ICP_DEFAULT_COS_NORMAL_THRESHOLD;
    min_distance_sq = (T)ICP_DEFAULT_MIN_DISTANCE_SQ;
    max_distance_sq = (T)ICP_DEFAULT_MAX_DISTANCE_SQ;
    icp_method = ICP_DEFAULT_METHOD;
  }

  template <typename T>
  void ICP<T>::match(Mat4x4<T>& ret_pc1_pc2, const T* pc1, 
    const uint32_t len_pc1, const T* pc2, const uint32_t len_pc2, 
    const Mat4x4<T>& guess_pc1_pc2, const T* norm_pc1, 
    const T* norm_pc2) {

    // Resize data structures if necessary
    if (transforms_.capacity() < num_iterations) {
      transforms_.capacity(num_iterations);
    }
    transforms_.resize(0);

    if (matches_.capacity() < len_pc2) {
      matches_.capacity(len_pc2);
    }
    matches_.resize(len_pc2);

    if (weights_.capacity() < len_pc2) {
      weights_.capacity(len_pc2);
    }
    weights_.resize(len_pc2);

    if (pc2_transformed_.capacity() / 3 < len_pc2) {
      pc2_transformed_.capacity(len_pc2 * 3);
    }
    pc2_transformed_.resize(len_pc2 * 3);

    if (norm_pc2 != NULL) {
      if (norm_pc2_transformed_.capacity() / 3 < len_pc2) {
        norm_pc2_transformed_.capacity(len_pc2 * 3);
      }
      norm_pc2_transformed_.resize(len_pc2 * 3);
    } else {
      norm_pc2_transformed_.resize(0);
    }
    
    transforms_.pushBack(guess_pc1_pc2);
    Mat4x4<T> cur_icp_mat;
    Mat4x4<T> new_composite_mat;
    for (uint32_t i = 0; i < num_iterations; i++) {
      if (verbose) {
        std::cout << "ICP iteration " << (i + 1) << " of " << num_iterations;
        std::cout << std::endl;
      }
      transformPC(&pc2_transformed_[0], &norm_pc2_transformed_[0], transforms_[i], 
        pc2, norm_pc2, len_pc2);
      if (norm_pc1 != NULL) {
        calcICPMat(cur_icp_mat, pc1, norm_pc1, len_pc1, &pc2_transformed_[0], 
          &norm_pc2_transformed_[0], len_pc2);
      } else {
        calcICPMat(cur_icp_mat, pc1, NULL, len_pc1, &pc2_transformed_[0], 
          NULL, len_pc2);
      }

      Mat4x4<T>::mult(new_composite_mat, cur_icp_mat, transforms_[i]);
      transforms_.pushBack(new_composite_mat);
    }
    ret_pc1_pc2.set(transforms_[transforms_.size() - 1]);
  }

  template <typename T>
  void ICP<T>::calcICPMat(Mat4x4<T>& ret, const T* D, const T* norm_D,
    const uint32_t len_D, const T* Q, const T* norm_Q, 
    const uint32_t len_Q) {
    // Note D is pc1 in the top level code, and Q is pc2

    // Zero out the matches and weights
    for (uint32_t i = 0; i < len_Q; i++) {
      matches_[i] = -1;
      weights_[i] = 1.0f;
    }

    // Construct a KD tree using pc2
    // Flann library expects features to be row-major (one point on each row)
    // luckily this is exactly how we have it.
    if (verbose) {
      std::cout << "    Finding correspondances..." << std::endl;
    }
    const int dim = 3;
    const int nn = 1;  // Num nearest neighbours to search for

    PointCloud<T> pc1(const_cast<T*>(D), len_D);
    PointCloud<T> pc2(const_cast<T*>(Q), len_Q);

    // construct a kd-tree index:
	  typedef KDTreeSingleIndexAdaptor<
		  L2_Simple_Adaptor<T, PointCloud<T> > ,
		  PointCloud<T>,
		  3 /* dim */
		  > my_kd_tree_t;

    my_kd_tree_t index(3 /*dim*/, pc1, 
      KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
	  index.buildIndex();

    const size_t num_results = 1;
    for (uint32_t i = 0; i < len_Q; i++) {
      size_t closest_index;
      T closest_distance_sq;
      index.knnSearch(&Q[i*3], num_results, &closest_index, 
        &closest_distance_sq);
      matches_[i] = static_cast<int>(closest_index);
      weights_[i] = static_cast<T>(closest_distance_sq);
    }

    /*
    // This is the old FLANN code.  I've since switched to nanoflann to make 
    // the code more portable (nanoflann is just a header library).
    flann::Matrix<T> dataset(const_cast<T*>(D), len_D, dim);
    flann::Matrix<T> query(const_cast<T*>(Q), len_Q, dim);

    flann::Matrix<int> indices(&matches_[0], len_Q, nn);
    flann::Matrix<T> dists(&weights_[0], len_Q, nn);

    // construct an randomized kd-tree index using 4 kd-trees
    // Note: L2 is the squared distances
    flann::Index<flann::L2<T>> index(dataset, flann::KDTreeIndexParams(4));
    index.buildIndex();

    // do a knn search, using 128 checks
    // For each point in the query, what is its closest neightbour in the dataset
    //flann::SearchParams sp = flann::SearchParams(128);
    flann::SearchParams sp = flann::SearchParams(32);
    sp.cores = ICP_NUM_CORES;  // Must be const in order for OpenMP to compile
                               // Multithreaded code!
    index.knnSearch(query, indices, dists, nn, sp);
    */

    // Compute the weights per match
    if (verbose) {
      std::cout << "    Calculating correspondance weights..." << std::endl;
    }

#if defined(DEBUG) || defined(_DEBUG)
    // Validate the FLANN min match
    uint32_t i_min = 0;
    T w_min = weights_[0];
    for (uint32_t i = 0; i < len_Q; i++) {
      if (weights_[i] < w_min) {
        w_min = weights_[i];
        i_min = i;
      }
    }
    Vec3<T> Q_min(&Q[i_min * 3]);
    Vec3<T> D_min(&D[matches_[i_min] * 3]);
    Vec3<T> diff;
    Vec3<T>::sub(diff, Q_min, D_min);
    T weight = Vec3<T>::dot(diff, diff);
    if (fabs(weight - w_min) > LOOSE_EPSILON) {
      std::cout << "WARNING! manually calculated weight doesn't match Flann's" << std::endl;
      std::cout << weight << " vs " << w_min << std::endl;
    }
    for (uint32_t i = 0; i < len_D; i++) {
      Vec3<T> D_val(&D[i * 3]);
      Vec3<T> diff;
      Vec3<T>::sub(diff, Q_min, D_val);
      weight = Vec3<T>::dot(diff, diff);
      if (weight < w_min) {
        std::cout << "WARNING! manually calculated weight is lower than Flann's" << std::endl;
        std::cout << weight << " vs " << w_min << std::endl;
      }
    }
#endif

    // weights_[i] is the squared distance of the correspondance
    T min_distance = sqrt(min_distance_sq);
    for (uint32_t i = 0; i < weights_.size(); i++) {
      if (weights_[i] >= max_distance_sq) {
        weights_[i] = 0.0f;
      } else {
        T angle_adjustment = 1.0f;
        if (norm_Q != NULL) {
          Vec3<T> cur_norm_D(&(norm_D[matches_[i] * 3]));
          Vec3<T> cur_norm_Q(&(norm_Q[i * 3]));
          T cos_angle = Vec3<T>::dot(cur_norm_D, cur_norm_Q);
          angle_adjustment = std::max<T>(0, 
            (cos_angle - cos_normal_threshold) / (1 - cos_normal_threshold));
        }
#ifdef ICP_LINEAR_WEIGHT_FUNCTION
        weights_[i] = std::max<T>(sqrt(weights_[i]), min_distance);
#else
        weights_[i] = std::max<T>(weights_[i], min_distance_sq);
#endif
        weights_[i] = angle_adjustment / (1.0f + weights_[i]);
      }
    }

    // Compute total weight and scale the weights
    T total_weight = 0.0f;
    for (uint32_t i = 0; i < weights_.size(); i++) {
      total_weight += weights_[i];
    }
    for (uint32_t i = 0; i < weights_.size(); i++) {
      weights_[i] /= total_weight;
    }

    switch (icp_method) {
    case SVD_ICP:
      {
        // Compute Mean of both sets
        if (verbose) {
          std::cout << "    Computing point cloud means..." << std::endl;
        }
        Vec3<T> D_mean, Q_mean;
        D_mean.zeros();
        Q_mean.zeros();
        uint32_t n_pts = 0;
        for (uint32_t i = 0; i < weights_.size(); i++) {
          if (weights_[i] > (T)EPSILON) {
            // Always true, but lets be explicit, just in case this changes
            n_pts++;
            Q_mean[0] += Q[i * 3] * weights_[i];
            Q_mean[1] += Q[i * 3 + 1] * weights_[i];
            Q_mean[2] += Q[i * 3 + 2] * weights_[i];
            uint32_t j = matches_[i];
            D_mean[0] += D[j * 3] * weights_[i];
            D_mean[1] += D[j * 3 + 1] * weights_[i];
            D_mean[2] += D[j * 3 + 2] * weights_[i];
          }
        }
        // No need to normalize since the weights have already been scaled

        // Calculate Cross-Covariance matrix
        if (verbose) {
          std::cout << "    Calculating Cross-Covariance Matrix..." << std::endl;
        }
        Vec3<T> D_pt, Q_pt;
        Mat3x3<T> cross_cov, outer_prod;
        for (uint32_t i = 0; i < weights_.size(); i++) {
          if (weights_[i] > (T)EPSILON) {
            Q_pt.set(&Q[i * 3]);
            D_pt.set(&D[matches_[i] * 3]);
            Vec3<T>::outerProd(outer_prod, Q_pt, D_pt);
            Mat3x3<T>::scale(outer_prod, weights_[i]);
            Mat3x3<T>::add(cross_cov, cross_cov, outer_prod);
          }
        }
        Vec3<T>::outerProd(outer_prod, Q_mean, D_mean);
        Mat3x3<T>::sub(cross_cov, cross_cov, outer_prod);

        // Compute SVD using Eigen
        if (verbose) {
          std::cout << "    Performing SVD..." << std::endl;
        }
        for (uint32_t i = 0; i < 3; i++) {
          for (uint32_t j = 0; j < 3; j++) {
            edata_->cross_cov_mat_(i, j) = (double)cross_cov(i, j);
          }
        }
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(edata_->cross_cov_mat_, 
          Eigen::ComputeFullU | Eigen::ComputeFullV);
        edata_->rot_e_mat_ = svd.matrixU() * svd.matrixV().transpose();
        if (edata_->rot_e_mat_.determinant() < 0) {
          Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
          D(2,2) = -1;
          edata_->rot_e_mat_ = svd.matrixU() * D * svd.matrixV().transpose();
          cout << "    fixing reflection" << endl;
        }
        Mat4x4<T> new_rot;
        new_rot.identity();
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            new_rot(i, j) = (T)edata_->rot_e_mat_(i,j);
          }
        }

        Mat4x4<T> T1;
        Mat4x4<T>::translationMat(T1, -Q_mean[0], -Q_mean[1], -Q_mean[2]);
        Mat4x4<T> T2;
        Mat4x4<T>::inverse(T2, new_rot);
        Mat4x4<T> T3;
        Mat4x4<T>::translationMat(T3, D_mean[0], D_mean[1], D_mean[2]);

        Mat4x4<T> temp;
        Mat4x4<T>::mult(temp, T2, T1);
        Mat4x4<T>::mult(ret, T3, temp);
      }
      break;
    case UMEYAMA_ICP:
      {
        if (verbose) {
          std::cout << "    Calculating Umeyama method matrix..." << std::endl;
        }
        // Perform Umeyama method
        uint32_t n_pts = 0;
        for (uint32_t i = 0; i < weights_.size(); i++) {
          if (weights_[i] > (T)EPSILON) {
            n_pts++;
          }
        }
        if (n_pts <= 3) {
          throw std::wruntime_error("ERROR: Not enough correspondance points!");
        }

        Eigen::MatrixXd X;  // Point set X is brought onto Y
        Eigen::MatrixXd Y;
        X.resize(3, n_pts);  // Each column is a point
        Y.resize(3, n_pts);

        // Fill up the Eigen structure
        uint32_t dst_ind = 0;
        for (uint32_t i = 0; i < matches_.size(); i++) {
          if (weights_[i] > (T)EPSILON) {
            X.col(dst_ind) <<  (double)Q[i * 3], (double)Q[i * 3 + 1], (double)Q[i * 3 + 2];
            uint32_t j = matches_[i];
            Y.col(dst_ind) << (double)D[j * 3], (double)D[j * 3 + 1], (double)D[j * 3 + 2];
            dst_ind++;
          }
        }
        Eigen::MatrixXd mat = Eigen::umeyama(X, Y, true);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            ret(i, j) = (T)mat(i,j);
          }
        }
        break;
      }
    case BFGS_ICP:
      {
        if (verbose) {
          std::cout << "    Calculating BFGS method matrix..." << std::endl;
        }

        // Collect the non-zero weighted correspondances and convert them to
        // double for numerical stability:
        uint32_t n_pts = 0;
        for (uint32_t i = 0; i < matches_.size(); i++) {
          if (weights_[i] > (T)EPSILON) {
            n_pts++;
          }
        }
        if (n_pts <= 3) {
          throw std::wruntime_error("ERROR: Not enough correspondance points ("
            "at least 3 is needed)!");
        }

        if (cur_Q_.capacity() < n_pts * 3) {
          cur_Q_.capacity(n_pts * 3);
        }
        cur_Q_.resize(n_pts * 3);
        if (cur_D_.capacity() < n_pts * 3) {
          cur_D_.capacity(n_pts * 3);
        }
        cur_D_.resize(n_pts * 3);
        if (cur_weights_.capacity() < n_pts) {
          cur_weights_.capacity(n_pts);
        }
        cur_weights_.resize(n_pts);
        uint32_t ind = 0;
        for (uint32_t i = 0; i < matches_.size(); i++) {
          if (weights_[i] > (T)EPSILON) {
            uint32_t j = matches_[i];
            cur_D_[ind * 3] = (double)D[j * 3];
            cur_D_[ind * 3 + 1] = (double)D[j * 3 + 1];
            cur_D_[ind * 3 + 2] = (double)D[j * 3 + 2];
            cur_Q_[ind * 3] = (double)Q[i * 3];
            cur_Q_[ind * 3 + 1] = (double)Q[i * 3 + 1];
            cur_Q_[ind * 3 + 2] = (double)Q[i * 3 + 2];
            cur_weights_[ind] = (double)weights_[i];
            ind++;
          }
        }

        // The starting coeffs should result in the identity matrix
        double start_coeff[ICPBFGSCoeffs::ICP_BFGS_NUM_COEFFS];
        for (uint32_t i = ICP_POS_X; i <= ICP_ORIENT_Z; i++) {
          start_coeff[i] = 0.0;
        }
#ifdef ICP_BFGS_INCLUDE_SCALE
        for (uint32_t i = ICP_SCALE_X; i <= ICP_SCALE_Z; i++) {
          start_coeff[i] = 2.0 * M_PI;
        }
#endif
        double end_coeff[ICPBFGSCoeffs::ICP_BFGS_NUM_COEFFS];
        bfgs_->verbose = verbose;
        bfgs_->c1 = 1e-8;
        bfgs_->descent_cond = ARMIJO;
        bfgs_->max_iterations = 100;
        bfgs_->jac_2norm_term = 1e-10;
        bfgs_->delta_f_term = 1e-12;
        bfgs_->delta_x_2norm_term = 1e-12;
        bfgs_->minimize(end_coeff, start_coeff, bfgs_angle_coeffs_, 
          bfgsObjFunc, bfgsJacobFunc, bfgsUpdateFunc);
        Double4x4 bfgs_ret;
        bfgsCoeffsToMat(bfgs_ret, end_coeff);
        for (uint32_t i = 0; i < 16; i++) {
          ret[i] = (T)bfgs_ret[i];
        }
        if (verbose) {
          std::cout << "    BFGS complete using pos, euler, and scale coeffs:";
          std::cout << std::endl;
          for (uint32_t i = 0; i < ICP_BFGS_NUM_COEFFS; i++) {
            std::cout << end_coeff[i] << ", ";
          }
          std::cout << std::endl;
        }
      }
      break;
    default:
      throw std::wruntime_error("ICP::match() - ERROR: ICPMethod is not "
        "recognized!");
    }
  }

  template <typename T>
  void ICP<T>::transformPC(T* pc_dst, T* norm_pc_dst, const Mat4x4<T>& mat, 
    const T* pc_src, const T* norm_pc_src, const uint32_t len_pc) {
    if (verbose) {
      std::cout << "    Transforming point cloud..." << std::endl;
    }
    Vec3<T> pt;
    Vec3<T> pt_transformed;
    for (uint32_t i = 0; i < len_pc; i++) {
      pt.set(&pc_src[i*3]);
      Vec3<T>::affineTransformPos(pt_transformed, mat, pt);
      pc_dst[i * 3] = pt_transformed[0];
      pc_dst[i * 3 + 1] = pt_transformed[1];
      pc_dst[i * 3 + 2] = pt_transformed[2];
    }

    // Now transform the normals
    if (norm_pc_src != NULL) {
      // In general ICP deals with only translation and rotation matricies,
      // so a full normal matrix is not really required.  However, since the
      // current matrix is seeded by an user-supplied transform we cannot 
      // assume that it is even affine, let alone only trans + rotation
      Mat4x4<T> normal_mat;
      Mat4x4<T>::inverse(normal_mat, mat);
      normal_mat.transpose();
      Vec3<T> norm;
      Vec3<T> norm_transformed;
      for (uint32_t i = 0; i < len_pc; i++) {
        norm.set(&norm_pc_src[i*3]);
        Vec3<T>::affineTransformVec(norm_transformed, normal_mat, norm);
        norm_transformed.normalize();
        norm_pc_dst[i * 3] = norm_transformed[0];
        norm_pc_dst[i * 3 + 1] = norm_transformed[1];
        norm_pc_dst[i * 3 + 2] = norm_transformed[2];
      }
    }
  }

  template <typename T>
  bool ICP<T>::bfgs_angle_coeffs_[ICP_BFGS_NUM_COEFFS] = {
    false,  // ICP_POSX
    false,  // ICP_POSY
    false,  // ICP_POSZ
    true,  // ICP_ORIENT_X
    true,  // ICP_ORIENT_Y
    true,  // ICP_ORIENT_Z
#ifdef ICP_BFGS_INCLUDE_SCALE
    false,  // ICP_SCALE_X
    false,  // ICP_SCALE_Y
    false,  // ICP_SCALE_Z
#endif
  };

  template <typename T>
  void ICP<T>::bfgsCoeffsToMat(Double4x4& mat, const double* coeff) {
    Double4x4::euler2RotMat(mat, coeff[ICP_ORIENT_X], 
      coeff[ICP_ORIENT_Y], coeff[ICP_ORIENT_Z]);
#ifdef ICP_BFGS_INCLUDE_SCALE
    mat.rightMultScale(coeff[ICP_SCALE_X] / (2.0 * M_PI), 
      coeff[ICP_SCALE_Y] / (2.0 * M_PI), coeff[ICP_SCALE_Z] / (2.0 * M_PI));
#endif
    mat.leftMultTranslation(coeff[ICP_POS_X] * 100.0, 
      coeff[ICP_POS_Y] * 100.0, coeff[ICP_POS_Z] * 100.0);
  }

  template <typename T>
  double ICP<T>::bfgsObjFunc(const double* coeff) {
    bfgsCoeffsToMat(cur_mat_, coeff);
    // Now calculate:
    // 1/N * Sum_i ( weight_i * ||D_i - mat * Q_i||_2 )
    
    // Now calculate the final objective function value
    Double3 D_pt, Q_pt, Q_pt_transformed, delta;
    double sum = 0.0;
    double scale = 1.0 / (((double)cur_D_.size()) / 3.0);
    for (uint32_t i = 0; i < cur_D_.size() / 3; i++) {
      D_pt.set(cur_D_[i * 3], cur_D_[i * 3 + 1], cur_D_[i * 3 + 2]);
      Q_pt.set(cur_Q_[i * 3], cur_Q_[i * 3 + 1], cur_Q_[i * 3 + 2]);
      Double3::affineTransformPos(Q_pt_transformed, cur_mat_, Q_pt);
      Double3::sub(delta, D_pt, Q_pt_transformed);
      sum += scale * cur_weights_[i] * Double3::dot(delta, delta);
    }

    // Add a penalty term for small scales so that the global solution of 
    // scale = 0, 0, 0, isn't a solution:

#ifdef ICP_BFGS_INCLUDE_SCALE
    sum += coeff[ICP_SCALE_X] <= 0.4 ? (1.0 / coeff[ICP_SCALE_X]) : 0.0;
    sum += coeff[ICP_SCALE_Y] <= 0.4 ? (1.0 / coeff[ICP_SCALE_Y]) : 0.0;
    sum += coeff[ICP_SCALE_Z] <= 0.4 ? (1.0 / coeff[ICP_SCALE_Z]) : 0.0;
#endif

    return sum;
  }

  template <typename T>
  void ICP<T>::bfgsJacobFunc(double* jacob, const double* coeff) {
    double coeff_tmp[ICP_BFGS_NUM_COEFFS];
    // Estimate using central diff.
    // http://math.fullerton.edu/mathews/n2003/differentiation/NumericalDiffProof.pdf
    memcpy(coeff_tmp, coeff, sizeof(coeff_tmp[0]) * ICP_BFGS_NUM_COEFFS);
    const double h = 0.00001;
    for (uint32_t i = 0; i < ICP_BFGS_NUM_COEFFS; i++) {
      coeff_tmp[i] = coeff[i] - h;
      double f0 = bfgsObjFunc(coeff_tmp);
      coeff_tmp[i] = coeff[i] + h;
      double f1 = bfgsObjFunc(coeff_tmp);
      coeff_tmp[i] = coeff[i];
      jacob[i] = (f1 - f0) / (2.0 * h);
    }
  }

  template <typename T>
  void ICP<T>::bfgsUpdateFunc(double* coeff) {
    WrapTwoPI(coeff[ICP_ORIENT_X]);
    WrapTwoPI(coeff[ICP_ORIENT_Y]);
    WrapTwoPI(coeff[ICP_ORIENT_Z]);
#ifdef ICP_BFGS_INCLUDE_SCALE
    coeff[ICP_SCALE_X] = fabs(coeff[ICP_SCALE_X]);
    coeff[ICP_SCALE_Y] = fabs(coeff[ICP_SCALE_Y]);
    coeff[ICP_SCALE_Z] = fabs(coeff[ICP_SCALE_Z]);
#endif
  }

  template <typename T>
  ICP<T>::~ICP() {
    SAFE_DELETE(edata_);
    SAFE_DELETE(bfgs_);
  }

}  // namespace math
}  // namespace icp

// Explicit template instantiation
template class icp::math::ICP<float>;
template class icp::math::ICP<double>;