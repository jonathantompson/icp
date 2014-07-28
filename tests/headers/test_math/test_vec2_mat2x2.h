//
//  test_vec2_mat2x2.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "test_unit/test_unit.h"
#include "icp/math/math_types.h"
#include "icp/math/math_base.h"

using icp::math::Vec2;
using icp::math::Mat2x2;

TEST(Vec2_Mat2x2, SimpleManipulation) {
  // Some practice variables
  Vec2<double> V_1(5.376671395461000e-001, 1.833885014595087e+000);
  Vec2<double> V_2(-2.258846861003648e+000, 8.621733203681206e-001);
#ifdef ROW_MAJOR
  Mat2x2<double> M_1(3.187652398589808e-001, -4.335920223056836e-001,
    -1.307688296305273e+000,  3.426244665386499e-001);
  Mat2x2<double> M_2(3.578396939725761e+000, -1.349886940156521e+000,
    2.769437029884877e+000,  3.034923466331855e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat2x2<double> M_1(3.187652398589808e-001, -1.307688296305273e+000, 
    -4.335920223056836e-001, 3.426244665386499e-001);
  Mat2x2<double> M_2(3.578396939725761e+000, 2.769437029884877e+000,  
    -1.349886940156521e+000, 3.034923466331855e+000);
#endif

  double M_2_DET = Mat2x2<double>::det(M_2);
  EXPECT_APPROX_EQ(M_2_DET, 1.459858772245127e+001);

  M_2.transpose();
#ifdef ROW_MAJOR
  Mat2x2<double> M_2_expect(3.578396939725761e+000, 2.769437029884877e+000,
    -1.349886940156521e+000, 3.034923466331855e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat2x2<double> M_2_expect(3.578396939725761e+000, -1.349886940156521e+000, 
    2.769437029884877e+000, 3.034923466331855e+000);
#endif
  EXPECT_TRUE(Mat2x2<double>::equal(M_2, M_2_expect));

  Mat2x2<double> M_3;
  Mat2x2<double>::mult(M_3, M_1, M_2);
#ifdef ROW_MAJOR
  Mat2x2<double> M_3_expect(1.725968767068822e+000, -4.331183442042076e-001,
    -5.141932090372603e+000, -2.581721357697310e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat2x2<double> M_3_expect(1.725968767068822e+000, -5.141932090372603e+000, 
    -4.331183442042076e-001, -2.581721357697310e+000);
#endif
  EXPECT_TRUE(Mat2x2<double>::equal(M_3, M_3_expect));

  Vec2<double> V_3;
  Vec2<double>::add(V_3, V_1, V_2);
  Vec2<double> V_3_expect(-1.721179721457548e+000, 2.696058334963207e+000);
  EXPECT_TRUE(Vec2<double>::equal(V_3, V_3_expect));

  Vec2<double> V_4;
  Vec2<double>::sub(V_4, V_1, V_2);
  Vec2<double> V_4_expect(2.796514000549748e+000, 9.717116942269659e-001);
  EXPECT_TRUE(Vec2<double>::equal(V_4, V_4_expect));

  Vec2<double> V_5;
  Vec2<double>::pairwiseMult(V_5, V_1, V_2);
  Vec2<double> V_5_expect(-1.214507730428518e+000, 1.581126732206785e+000);
  EXPECT_TRUE(Vec2<double>::equal(V_5, V_5_expect));

  double VOT = Vec2<double>::dot(V_1, V_2);
  EXPECT_APPROX_EQ(VOT, 3.666190017782667e-001);

  Mat2x2<double> M_4;
  Mat2x2<double>::inverse(M_4, M_3);
#ifdef ROW_MAJOR
  Mat2x2<double> M_4_expect(3.863096854126326e-001, -6.480862498856263e-002, 
    -7.694006800240641e-001, -2.582611982700881e-001);
#endif
#ifdef COLUMN_MAJOR
  Mat2x2<double> M_4_expect(3.863096854126326e-001, -7.694006800240641e-001, 
    -6.480862498856263e-002, -2.582611982700881e-001);
#endif
  EXPECT_TRUE(Mat2x2<double>::equal(M_4, M_4_expect));
}
