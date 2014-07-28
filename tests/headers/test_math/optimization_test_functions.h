//
//  optimization_test_functions.h
//
//  Created by Jonathan Tompson on 10/15/12.
//
//  Some functions for testing optimization routines.
//
//  You must call 'initOptimizationTestFunctions()' before using these functions
//

#ifndef ik_MATH_OPTIMIZATION_TEST_FUNCTIONS_HEADER
#define ik_MATH_OPTIMIZATION_TEST_FUNCTIONS_HEADER

#include "icp/math/math_types.h"

#define NUM_COEFFS_ROSENBROCK 4  // Don't change!
#define NUM_COEFFS_RASTRIGIN 4  // Don't change!
#define NUM_COEFFS_ROTELLIPS 4  // Don't change!
#define NUM_COEFFS_EXPONTIAL_FIT 4  // Don't change!
#define NUM_PTS_EXPONTIAL_FIT 100  // Don't change!
#define NUM_COEFFS_HW7_4A 3  // Don't change!
#define NUM_COEFFS_HW7_4B 4  // Don't change!

#define C_DIM_HW3_3 4  // Don't change!
#define X_DIM_HW3_3 1  // Don't change!
#define NUM_PTS_HW3_3 100  // Don't change!

namespace icp {
namespace math {

  // Extended Rosenbrock function
  extern float c_0_rosenbrock[NUM_COEFFS_ROSENBROCK];
  extern float c_answer_rosenbrock[NUM_COEFFS_ROSENBROCK];
  float extendedRosenbrock(const float* coeff);

  // Generalized Rastragin function
  extern float c_0_rastrigin[NUM_COEFFS_ROSENBROCK];
  extern float c_answer_rastrigin[NUM_COEFFS_ROSENBROCK];
  float generalizedRastrigin(const float* coeff);

  // Rotated ellipsoidal function
  extern float c_0_rotellips[NUM_COEFFS_ROTELLIPS];
  extern float c_answer_rotellips[NUM_COEFFS_ROTELLIPS];
  float rotatedEllipsoidal(const float* coeff);

  // Exponential fit function (see matlab code non_linear_least_sq.m)
  extern float y_vals_rand[NUM_PTS_EXPONTIAL_FIT];
  extern float x_vals_rand[NUM_PTS_EXPONTIAL_FIT];
  extern float f_exponential_fit[NUM_PTS_EXPONTIAL_FIT];
  extern float c_0_exponential_fit[NUM_COEFFS_EXPONTIAL_FIT];
  extern float c_answer_exponential_fit[NUM_COEFFS_EXPONTIAL_FIT];
  float exponentialFit(const float* coeff);

  // Non-linear R^3 function from HW7 (Q4a) of Numerical Optimization - NYU
  extern float c_0_hw7_4a[NUM_COEFFS_HW7_4A];
  extern float c_answer_hw7_4a[NUM_COEFFS_HW7_4A];
  float hw7_4a(const float* coeff);
  void hw7_4a_jacob(float* jacob, const float* coeff);
  
  extern float c_0_hw7_4b[NUM_COEFFS_HW7_4B];
  extern float c_answer_hw7_4b[NUM_COEFFS_HW7_4B];
  float hw7_4b(const float* coeff);
  void hw7_4b_jacob(float* jacob, const float* coeff);

  extern double dc_0_hw7_4b[NUM_COEFFS_HW7_4B];
  extern double dc_answer_hw7_4b[NUM_COEFFS_HW7_4B];
  double dhw7_4b(const double* coeff);
  void dhw7_4b_jacob(double* jacob, const double* coeff);

  extern double dc_start_hw3_3[C_DIM_HW3_3];
  extern double dc_answer_hw3_3[C_DIM_HW3_3];
  extern double dy_vals_hw3_3[NUM_PTS_HW3_3];
  extern double dx_vals_hw_3_3[NUM_PTS_HW3_3];
  double dfunc_hw3_3(const double* x, const double* c);
  void djacob_hw3_3(double* jacob, const double* x, const double* c);

  extern float c_start_hw3_3[C_DIM_HW3_3];
  extern float c_answer_hw3_3[C_DIM_HW3_3];
  extern float y_vals_hw3_3[NUM_PTS_HW3_3];
  extern float x_vals_hw_3_3[NUM_PTS_HW3_3];
  float func_hw3_3(const float* x, const float* c);
  void jacob_hw3_3(float* jacob, const float* x, const float* c);

}  // namespace math
}  // namespace icp

#endif  // ik_MATH_OPTIMIZATION_TEST_FUNCTIONS_HEADER
