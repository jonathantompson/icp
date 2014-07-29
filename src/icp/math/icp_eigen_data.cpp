#include <random>
#include <stdexcept>
#include <iostream>
#include "icp/math/icp_eigen_data.h"

using Eigen::MatrixXd;
using std::cout;
using std::endl;
using std::runtime_error;

namespace icp {
namespace math {
  ICPEigenData::ICPEigenData() {
    cross_cov_mat_.resize(3, 3);
    rot_e_mat_.resize(3, 3);
  }

  ICPEigenData::~ICPEigenData() {

  }

}  // namespace math
}  // namespace icp
