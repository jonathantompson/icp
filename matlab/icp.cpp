#include <sstream>
#include "mex.h"
#include "matrix.h"
#include "icp/icp_base.h"

using namespace icp::math;

double* getFieldSafe(const mxArray* arr, const int index, const char* field_name) {
  mxArray* internal_array = mxGetField(arr, index, field_name);
  if (internal_array == NULL) {
    std::stringstream ss;
    ss << field_name << " field is missing";
    mexErrMsgIdAndTxt("MATLAB:icp:invalidInput", ss.str().c_str());
  }
  if (mxGetClassID(internal_array) != mxDOUBLE_CLASS) {
    std::stringstream ss;
    ss << field_name << " field is not double type";
    mexErrMsgIdAndTxt("MATLAB:icp:invalidInput", ss.str().c_str());
  }
  if (mxGetNumberOfDimensions(internal_array) != 2) {
    std::stringstream ss;
    ss << field_name << " has an invalid number of dimensions";
    mexErrMsgIdAndTxt("MATLAB:icp:invalidInput", ss.str().c_str());
  }
  return static_cast<double*>(mxGetData(internal_array));
}

double* getOptionalFieldSafe(const mxArray* arr, const int index, const char* field_name) {
  mxArray* internal_array = mxGetField(arr, index, field_name);
  if (internal_array == NULL || mxGetClassID(internal_array) != mxDOUBLE_CLASS) {
    return NULL;
  }
  return static_cast<double*>(mxGetData(internal_array));
}

template <class T>
T* getOptionalFieldSafe(const mxArray* arr, const int index, 
  const char* field_name, const mxClassID classId) {
  mxArray* internal_array = mxGetField(arr, index, field_name);
  if (internal_array == NULL) {
    return NULL;
  }
  if (mxGetClassID(internal_array) != classId) {
    return NULL;
  }
  return (T*)mxGetData(internal_array);
}

// The gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // First check that the inputs and output exist
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:icp:invalidNumOutput",
      "You must define one (and only one) return value!");
  }
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:icp:invalidNumInput",
      "Invalid number of inputs");
  } 
  if (!mxIsStruct(prhs[0])) {
    mexErrMsgIdAndTxt("MATLAB:icp:invalidInput",
      "Invalid input type (struct expected)");
  }   
  
  // Get the fields
  double* p_num_iterations = getFieldSafe(prhs[0], 0, "num_iterations");
  double* p_pc1 = getFieldSafe(prhs[0], 0, "pc1");
  double* p_m_pc2_initial = getFieldSafe(prhs[0], 0, "m_pc2_initial");
  double* p_pc2 = getFieldSafe(prhs[0], 0, "pc2");
  double* p_min_distance_sq = getOptionalFieldSafe(prhs[0], 0, "min_distance_sq");
  double* p_max_distance_sq = getOptionalFieldSafe(prhs[0], 0, "max_distance_sq");
  double* p_cos_normal_threshold = getOptionalFieldSafe(prhs[0], 0, 
    "cos_normal_threshold");
  double* p_npc1 = getOptionalFieldSafe(prhs[0], 0, "npc1");
  double* p_npc2 = getOptionalFieldSafe(prhs[0], 0, "npc2");
  if ((p_npc1 == NULL) != (p_npc2 == NULL)) {
    mexErrMsgIdAndTxt("MATLAB:icp:invalidInput",
      "either define npc1 AND npc2 or don't define both");
  }
  double* p_method = getOptionalFieldSafe(prhs[0], 0, "method");
  double* p_match_scale = getOptionalFieldSafe(prhs[0], 0, "match_scale");
  
  // Check the dimension and sizing of pc1, pc2
  mxArray* arr;
  arr = mxGetField(prhs[0], 0, "pc1");
  uint32_t num_points_pc1 = static_cast<uint32_t>(mxGetN(arr));
  uint32_t dim_pc1 = static_cast<uint32_t>(mxGetM(arr));
  arr = mxGetField(prhs[0], 0, "pc2");
  uint32_t num_points_pc2 = static_cast<uint32_t>(mxGetN(arr));
  uint32_t dim_pc2 = static_cast<uint32_t>(mxGetM(arr));
  
  if (dim_pc1 != 3 || dim_pc2 != 3) {
    mexErrMsgIdAndTxt("MATLAB:icp:invalidInput",
      "pc1 or pc2 dimension is not 3 (number of rows should be 3)");
  }
  
  // Check the dimension and sizing of npc1, npc2
  if (p_npc1 != NULL) {
    arr = mxGetField(prhs[0], 0, "npc1");
    uint32_t num_points = static_cast<uint32_t>(mxGetN(arr));
    uint32_t dim = static_cast<uint32_t>(mxGetM(arr));
    if (num_points != num_points_pc1 || dim != 3) {
      mexErrMsgIdAndTxt("MATLAB:icp:invalidInput",
        "npc1 is not the correct size!");
    }
  }
  if (p_npc2 != NULL) {
    arr = mxGetField(prhs[0], 0, "npc2");
    uint32_t num_points = static_cast<uint32_t>(mxGetN(arr));
    uint32_t dim = static_cast<uint32_t>(mxGetM(arr));
    if (num_points != num_points_pc2 || dim != 3) {
      mexErrMsgIdAndTxt("MATLAB:icp:invalidInput",
        "npc2 is not the correct size!");
    }
  }
  
  ICP<double>* p_icp = new ICP<double>();
  p_icp->num_iterations = static_cast<uint32_t>(*p_num_iterations);
  if (p_min_distance_sq != NULL) {
    p_icp->min_distance_sq = static_cast<float>(*p_min_distance_sq);
  }
  if (p_max_distance_sq != NULL) {
    p_icp->max_distance_sq = static_cast<float>(*p_max_distance_sq);
  }
  if (p_cos_normal_threshold != NULL) {
    p_icp->cos_normal_threshold = static_cast<float>(*p_cos_normal_threshold);
  }
  p_icp->verbose = false;
  int method = (int)ICPMethod::BFGS_ICP;
  if (p_method != NULL) {
    method = static_cast<int>(*p_method);
  }
  if (method < 0 || method >= (int)ICPMethod::NUM_METHODS) {
    mexErrMsgIdAndTxt("MATLAB:icp:invalidInput",
      "method options are 0 - SVD, 1 - BFGS, 2 - PSO, 3 - Umeyama!");
  }
  p_icp->icp_method = (icp::math::ICPMethod)method;
  p_icp->match_scale = true;
  if (p_match_scale != NULL) {
	p_icp->match_scale = *p_match_scale != 0;
  }

  Mat4x4<double> m_pc2_initial;
  Mat4x4<double> m_pc2_final;
  memcpy(&m_pc2_final[0], p_m_pc2_initial, sizeof(m_pc2_final[0]) * 4 * 4);
  
  p_icp->match(m_pc2_final, p_pc1, num_points_pc1, p_pc2, num_points_pc2, 
    m_pc2_initial, p_npc1, p_npc2);

  delete p_icp;
  
  // Allocate the output array and copy the final matrix to the output
  plhs[0] = mxCreateDoubleMatrix(4, 4, mxREAL);
  double* m_pc2_return = mxGetPr(plhs[0]);
  for (uint32_t i = 0; i < 16; i++) {
    m_pc2_return[i] = m_pc2_final[i];
  }
}