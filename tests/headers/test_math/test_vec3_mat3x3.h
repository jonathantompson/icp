//
//  test_vec2_mat2x2.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "test_unit/test_unit.h"
#include "icp/math/math_types.h"
#include "icp/math/math_base.h"

using icp::math::Vec3;
using icp::math::Mat3x3;

TEST(Vec3_Mat3x3, SimpleManipulation) {
  // Some practice variables
  Vec3<double> V_1(5.376671395461000e-001, 1.833885014595087e+000,
    -2.258846861003648e+000);
  Vec3<double> V_2(8.621733203681206e-001, 3.187652398589808e-001,
    -1.307688296305273e+000);
#ifdef ROW_MAJOR
  Mat3x3<double> M_1(-4.335920223056836e-001, 2.769437029884877e+000, 
    7.254042249461056e-001, 3.426244665386499e-001, -1.349886940156521e+000, 
    -6.305487318965619e-002, 3.578396939725761e+000, 3.034923466331855e+000,
    7.147429038260958e-001);
  Mat3x3<double> M_2(-2.049660582997746e-001, 1.409034489800479e+000, 
    -1.207486922685038e+000, -1.241443482163119e-001, 1.417192413429614e+000, 
    7.172386513288385e-001, 1.489697607785465e+000, 6.714971336080805e-001, 
    1.630235289164729e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat3x3<double> M_1(-4.335920223056836e-001,  3.426244665386499e-001,
    3.578396939725761e+000, 2.769437029884877e+000, -1.349886940156521e+000,
    3.034923466331855e+000, 7.254042249461056e-001, -6.305487318965619e-002,
    7.147429038260958e-001);
  Mat3x3<double> M_2(-2.049660582997746e-001, -1.241443482163119e-001,
    1.489697607785465e+000, 1.409034489800479e+000, 1.417192413429614e+000,
    6.714971336080805e-001, -1.207486922685038e+000, 7.172386513288385e-001,
    1.630235289164729e+000);
#endif

  double M_2_DET = Mat3x3<double>::det(M_2);
  EXPECT_APPROX_EQ(M_2_DET, 4.065739015327806e+000);

  M_2.transpose();
#ifdef ROW_MAJOR
  Mat3x3<double> M_2_expect(-2.049660582997746e-001, -1.241443482163119e-001, 
    1.489697607785465e+000, 1.409034489800479e+000, 1.417192413429614e+000, 
    6.714971336080805e-001, -1.207486922685038e+000, 7.172386513288385e-001, 
    1.630235289164729e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat3x3<double> M_2_expect(-2.049660582997746e-001, 1.409034489800479e+000,
    -1.207486922685038e+000, -1.241443482163119e-001, 1.417192413429614e+000,
    7.172386513288385e-001, 1.489697607785465e+000, 6.714971336080805e-001,
    1.630235289164729e+000);
#endif
  EXPECT_TRUE(Mat3x3<double>::equal(M_2, M_2_expect));

  Mat3x3<double> M_3;
  Mat3x3<double>::mult(M_3, M_1, M_2);
#ifdef ROW_MAJOR
  Mat3x3<double> M_3_expect(3.115187824877718e+000, 4.498941095193409e+000, 
    2.396327595308463e+000, -1.896125707607191e+000, -2.000809813865199e+000, 
    -4.988326422662731e-001, 2.679819212746901e+000, 4.369473992470597e+000, 
    8.533870873675408e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat3x3<double> M_3_expect(3.115187824877718e+000, -1.896125707607191e+000, 
    2.679819212746901e+000, 4.498941095193409e+000, -2.000809813865199e+000, 
    4.369473992470597e+000, 2.396327595308463e+000, -4.988326422662731e-001,
    8.533870873675408e+000);
#endif
  EXPECT_TRUE(Mat3x3<double>::equal(M_3, M_3_expect));

  Vec3<double> V_3;
  Vec3<double>::add(V_3, V_1, V_2);
  Vec3<double> V_3_expect(1.399840459914221e+000, 2.152650254454067e+000, 
    -3.566535157308922e+000);
  EXPECT_TRUE(Vec3<double>::equal(V_3, V_3_expect));

  Vec3<double> V_4;
  Vec3<double>::sub(V_4, V_1, V_2);
  Vec3<double> V_4_expect(-3.245061808220205e-001, 1.515119774736106e+000, 
    -9.511585646983747e-001);
  EXPECT_TRUE(Vec3<double>::equal(V_4, V_4_expect));

  Vec3<double> V_5;
  Vec3<double>::pairwiseMult(V_5, V_1, V_2);
  Vec3<double> V_5_expect(4.635622629552907e-001, 5.845787965511933e-001, 
    2.953867603280375e+000);
  EXPECT_TRUE(Vec3<double>::equal(V_5, V_5_expect));

  double VOT = Vec3<double>::dot(V_1, V_2);
  EXPECT_APPROX_EQ(VOT, 4.002008662786859e+000);

  Mat3x3<double> M_4;
  Mat3x3<double>::inverse(M_4, M_3);
#ifdef ROW_MAJOR
  Mat3x3<double> M_4_expect(-1.113337365293549e+000, -2.087099122677174e+000, 
    1.906295402620305e-001, 1.109562285946668e+000, 1.507087590420958e+000,   
    -2.234730602210047e-001, -2.185011605272204e-001, -1.162581104809674e-001, 
    1.717400042533082e-001);
#endif
#ifdef COLUMN_MAJOR
  Mat3x3<double> M_4_expect(-1.113337365293549e+000, 1.109562285946668e+000,
    -2.185011605272204e-001, -2.087099122677174e+000, 1.507087590420958e+000,
    -1.162581104809674e-001, 1.906295402620305e-001, -2.234730602210047e-001,
    1.717400042533082e-001);
#endif
  EXPECT_TRUE(Mat3x3<double>::equal(M_4, M_4_expect));
}

TEST(Vec3_Mat3x3, SymMatEig) {
  // M_1 is symmetric and so is the same row-major or column-major 
  Mat3x3<double> M_1(1.220430829109408e+000, 1.112291359983057e+000, 
    -3.893526056568728e+000, 1.112291359983057e+000, 3.582137049969669e+000, 
    -3.333264639555196e+000, -3.893526056568728e+000, -3.333264639555196e+000, 
    1.961736247999851e+001);

  Vec3<double> eigVecs[3];
  Vec3<double> eigVals;
  icp::math::SymMatEigenVecs(eigVecs, &eigVals, &M_1);  // O(n^3)
  Vec3<double> eigVals_expect(3.605466872560469e-001, 2.957209170251527e+000,
    2.110217450157001e+001);
  Vec3<double> eigVecs_expect[3];
  eigVecs_expect[0].set(9.723990622249192e-001, -1.611732403643351e-001, 
    1.687105520545704e-001);
  eigVecs_expect[1].set(1.218240883706243e-001, 9.674011029452027e-001, 
    2.220225157795305e-001);
  eigVecs_expect[2].set(-1.989948624381140e-001, -1.953414769342800e-001, 
    9.603347083763922e-001);
  // Matlab may return eigenvalues out of order, it may also return eigen 
  // vectors that have been flipped (since a negative eigenvector is still an
  // eigenvector
  bool eig_correct;
  bool eig_vec_correct;
  for (int i = 0; i < 3; i ++) {
    eig_correct = false;
    eig_vec_correct = false;
    for (int j = 0; j < 3; j ++) {
      if (Vec3<double>::equal(eigVecs[i], eigVecs_expect[j]))
        eig_vec_correct = true;
      icp::math::Vec3<double>::scale(eigVecs[i], -1);
      if (Vec3<double>::equal(eigVecs[i], eigVecs_expect[j]))
        eig_vec_correct = true;
      if (abs(eigVals.m[i] - eigVals_expect[j]) < EPSILON)
        eig_correct = true;
    }
    if (eig_correct == false || eig_vec_correct == false)
      break;
  }
  EXPECT_TRUE(eig_correct);
  EXPECT_TRUE(eig_vec_correct);
  icp::math::SymMatEigenVals(&eigVals, &M_1);  // O(n^2)
  // Matlab may return eigenvalues out of order
  for (int i = 0; i < 3; i ++) {
    eig_correct = false;
    for (int j = 0; j < 3; j ++) {
      if (abs(eigVals.m[i] - eigVals_expect[j]) < EPSILON)
        eig_correct = true;
    }
    if (eig_correct == false)
      break;
  }
  EXPECT_TRUE(eig_correct);   
}
