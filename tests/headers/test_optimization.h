//
//  test_optimization.cpp
//
//  Created by Jonathan Tompson on 5/25/13.
//
//  Test BFGS, and other optimizers.
//

#include "icp/math/pso.h"
#include "icp/math/pso_parallel.h"
#include "icp/math/bfgs.h"
#include "icp/math/lm_fit.h"
#include "icp/data_str/vector.h"
#include "test_math/optimization_test_functions.h"

TEST(PSO, ExponentialFit) {
  icp::math::PSO<float>* solver = 
    new icp::math::PSO<float>(NUM_COEFFS_EXPONTIAL_FIT, 17);
  solver->max_iterations = 10000;
  solver->delta_coeff_termination = 1e-8f;

  float ret_coeffs[NUM_COEFFS_EXPONTIAL_FIT];
  float c_rad[NUM_COEFFS_EXPONTIAL_FIT] = {2, 2, 2, 2};

  solver->minimize(ret_coeffs, icp::math::c_0_exponential_fit, 
    c_rad, NULL, icp::math::exponentialFit, NULL);

  for (uint32_t i = 0; i < NUM_COEFFS_EXPONTIAL_FIT; i++) {
    float k = fabsf(ret_coeffs[i] - icp::math::c_answer_exponential_fit[i]);
    EXPECT_TRUE(k < 0.001f);
  }

  delete solver;
}

namespace icp {
namespace math {
  void coeffUpdateFunc(float* coeff) { 

  }

  void exponentialFitParallel(icp::data_str::Vector<float>& residues, 
    icp::data_str::Vector<float*>& coeffs) {
    for (uint32_t i = 0; i < coeffs.size(); i++) {
      residues[i] = exponentialFit(coeffs[i]);
    }
  }
}  // namespace math
}  // namespace icp

TEST(PSOParallel, ExponentialFit) {
  bool angle_coeffs[NUM_COEFFS_EXPONTIAL_FIT];
  memset(angle_coeffs, 0, sizeof(angle_coeffs[0]) * NUM_COEFFS_EXPONTIAL_FIT);
  icp::math::PSOParallel* solver2 = 
    new icp::math::PSOParallel(NUM_COEFFS_EXPONTIAL_FIT, 64);

  float ret_coeffs2[NUM_COEFFS_EXPONTIAL_FIT];
  float c_rad[NUM_COEFFS_EXPONTIAL_FIT] = {2, 2, 2, 2};

  solver2->minimize(ret_coeffs2, icp::math::c_0_exponential_fit, c_rad, 
    angle_coeffs, icp::math::exponentialFitParallel, 
    &icp::math::coeffUpdateFunc);

  for (uint32_t i = 0; i < NUM_COEFFS_EXPONTIAL_FIT; i++) {
    float k = fabsf(ret_coeffs2[i] - icp::math::c_answer_exponential_fit[i]);
    EXPECT_TRUE(k < 0.001f);
  }

  delete solver2;
}

TEST(BFGS, A_FLOAT) {
  icp::math::BFGS<float>* solver_bfgs = new icp::math::BFGS<float>(NUM_COEFFS_HW7_4A);
  solver_bfgs->verbose = false;
  solver_bfgs->max_iterations = 1000;
  solver_bfgs->delta_f_term = 1e-12f;
  solver_bfgs->jac_2norm_term = 1e-12f;
  solver_bfgs->delta_x_2norm_term = 1e-12f;

  float ret_coeffs_bfgs[NUM_COEFFS_HW7_4A];

  solver_bfgs->minimize(ret_coeffs_bfgs, icp::math::c_0_hw7_4a, NULL, 
    icp::math::hw7_4a, icp::math::hw7_4a_jacob, NULL);

  for (uint32_t i = 0; i < NUM_COEFFS_HW7_4A; i++) {
    float k = fabsf(ret_coeffs_bfgs[i] - icp::math::c_answer_hw7_4a[i]);
    EXPECT_TRUE(k < 0.00001f);
  }

  float f_bfgs = icp::math::hw7_4a(ret_coeffs_bfgs);
  float f_answer = icp::math::hw7_4a(icp::math::c_answer_hw7_4a);
  EXPECT_TRUE(fabsf(f_bfgs - f_answer) < 0.00001f);

  delete solver_bfgs;
}

TEST(BFGS, B_FLOAT) {
  icp::math::BFGS<float>* solver_bfgs2 = new icp::math::BFGS<float>(NUM_COEFFS_HW7_4B);
  solver_bfgs2->verbose = false;
  solver_bfgs2->descent_cond = icp::math::SufficientDescentCondition::ARMIJO;
  solver_bfgs2->max_iterations = 1000;
  solver_bfgs2->delta_f_term = 1e-12f;
  solver_bfgs2->jac_2norm_term = 1e-12f;
  solver_bfgs2->delta_x_2norm_term = 1e-12f;
  
  float ret_coeffs_bfgs2[NUM_COEFFS_HW7_4B];
  
  solver_bfgs2->minimize(ret_coeffs_bfgs2, icp::math::c_0_hw7_4b, NULL, 
    icp::math::hw7_4b, icp::math::hw7_4b_jacob, NULL);
  
  for (uint32_t i = 0; i < NUM_COEFFS_HW7_4B; i++) {
    float k = fabsf(ret_coeffs_bfgs2[i] - icp::math::c_answer_hw7_4b[i]);
    EXPECT_TRUE(k < 0.00001f);
  }

  float f_bfgs = icp::math::hw7_4b(ret_coeffs_bfgs2);
  float f_answer = icp::math::hw7_4b(icp::math::c_answer_hw7_4b);
  EXPECT_TRUE(fabsf(f_bfgs - f_answer) < 0.00001f);

  delete solver_bfgs2;
}

TEST(BFGS, B_DOUBLE) {
  icp::math::BFGS<double>* solver_bfgs2 = new icp::math::BFGS<double>(NUM_COEFFS_HW7_4B);
  solver_bfgs2->verbose = false;
  solver_bfgs2->descent_cond = icp::math::SufficientDescentCondition::STRONG_WOLFE;
  solver_bfgs2->max_iterations = 1000;
  solver_bfgs2->delta_f_term = 1e-12;
  solver_bfgs2->jac_2norm_term = 1e-12;
  solver_bfgs2->delta_x_2norm_term = 1e-12;
  
  double ret_coeffs_bfgs2[NUM_COEFFS_HW7_4B];
  
  solver_bfgs2->minimize(ret_coeffs_bfgs2, icp::math::dc_0_hw7_4b, NULL, 
    icp::math::dhw7_4b, icp::math::dhw7_4b_jacob, NULL);
  
  for (uint32_t i = 0; i < NUM_COEFFS_HW7_4B; i++) {
    double k = fabs(ret_coeffs_bfgs2[i] - icp::math::dc_answer_hw7_4b[i]);
    EXPECT_TRUE(k < 0.00001);
  }

  double f_bfgs = icp::math::dhw7_4b(ret_coeffs_bfgs2);
  double f_answer = icp::math::dhw7_4b(icp::math::dc_answer_hw7_4b);
  EXPECT_TRUE(fabs(f_bfgs - f_answer) < 0.00001);

  delete solver_bfgs2;
}

TEST(LM_FIT, DOUBLE) {
  icp::math::LMFit<double>* lm_fit = 
    new icp::math::LMFit<double>(C_DIM_HW3_3, X_DIM_HW3_3, NUM_PTS_HW3_3);
  lm_fit->verbose = false;
  lm_fit->delta_c_termination = 1e-16;
  
  double ret_coeffs[C_DIM_HW3_3];
  
  lm_fit->fitModel(ret_coeffs, icp::math::dc_start_hw3_3, 
    icp::math::dy_vals_hw3_3, icp::math::dx_vals_hw_3_3, 
    icp::math::dfunc_hw3_3, icp::math::djacob_hw3_3, NULL, NULL);
  
  for (uint32_t i = 0; i < C_DIM_HW3_3; i++) {
    double k = fabs(ret_coeffs[i] - icp::math::dc_answer_hw3_3[i]);
    EXPECT_TRUE(k < 0.00001);
  }

  delete lm_fit;
}

TEST(LM_FIT, FLOAT) {
  icp::math::LMFit<float>* lm_fit = 
    new icp::math::LMFit<float>(C_DIM_HW3_3, X_DIM_HW3_3, NUM_PTS_HW3_3);
  lm_fit->verbose = false;
  lm_fit->delta_c_termination = 1e-16f;
  
  float ret_coeffs[C_DIM_HW3_3];
  
  lm_fit->fitModel(ret_coeffs, icp::math::c_start_hw3_3, 
    icp::math::y_vals_hw3_3, icp::math::x_vals_hw_3_3, 
    icp::math::func_hw3_3, icp::math::jacob_hw3_3, NULL, NULL);
  
  for (uint32_t i = 0; i < C_DIM_HW3_3; i++) {
    double k = fabs(ret_coeffs[i] - icp::math::c_answer_hw3_3[i]); 
    EXPECT_TRUE(k < 0.00001);
  }

  delete lm_fit;
}