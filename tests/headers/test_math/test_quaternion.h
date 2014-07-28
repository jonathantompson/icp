//
//  test_vec2_mat2x2.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "test_unit/test_unit.h"
#include "icp/math/math_types.h"
#include "icp/math/math_base.h"

using icp::math::Quat;
using icp::math::Mat3x3;
using icp::math::Mat4x4;

TEST(Quaternion, SimpleManipulation) {
  Vec3<double> axis(1.817168119886437e-001, 6.198030601001485e-001, 
    -7.634285604633715e-001);
  double angle = 8.621733203681206e-001;

  Quat<double> Q1(axis, angle);  // axis angle -> quaternion
  Quat<double> Q1_expect(7.593187686640246e-002, 2.589898486876663e-001, 
    -3.190049551002597e-001, 9.085122161940182e-001);
  EXPECT_TRUE(Quat<double>::equal(Q1, Q1_expect));

  Mat3x3<double> RotMat1;
  Quat<double>::quat2Mat3x3(RotMat1, Q1);  // quaternion -> matrix
#ifdef ROW_MAJOR
  Mat3x3<double> RotMat1_expect(6.623201937964420e-001, 
    6.189709680704206e-001, 4.221455928650800e-001, -5.403086268696203e-001, 
    7.849403773940535e-001, -3.032081655673973e-001, -5.190361727468607e-001,
    -2.726801464073089e-002, 8.543172167045705e-001);
#endif
#ifdef COLUMN_MAJOR
  Mat3x3<double> RotMat1_expect(6.623201937964420e-001, 
    -5.403086268696203e-001, -5.190361727468607e-001, 6.189709680704206e-001,
    7.849403773940535e-001, -2.726801464073089e-002, 4.221455928650800e-001,
    -3.032081655673973e-001, 8.543172167045705e-001);
#endif
  EXPECT_TRUE(Mat3x3<double>::equal(RotMat1, RotMat1_expect));

  Mat4x4<double> RotMat2;
  Quat<double>::quat2Mat4x4(RotMat2, Q1);  // quaternion -> matrix
#ifdef ROW_MAJOR
  Mat4x4<double> RotMat2_expect(6.623201937964420e-001, 
    6.189709680704206e-001, 4.221455928650800e-001, 0, 
    -5.403086268696203e-001, 7.849403773940535e-001, -3.032081655673973e-001,
    0, -5.190361727468607e-001, -2.726801464073089e-002, 
    8.543172167045705e-001, 0, 0, 0, 0, 1);
#endif
#ifdef COLUMN_MAJOR
  Mat4x4<double> RotMat2_expect(6.623201937964420e-001, 
    -5.403086268696203e-001, -5.190361727468607e-001, 0, 
    6.189709680704206e-001, 7.849403773940535e-001, -2.726801464073089e-002, 
    0, 4.221455928650800e-001, -3.032081655673973e-001, 
    8.543172167045705e-001, 0, 0, 0, 0, 1);
#endif
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat2, RotMat2_expect));

  Quat<double> Q2(RotMat1);  // matrix --> quaternion
  Quat<double> Q2_expect(7.593187686640246e-002, 2.589898486876663e-001, 
    -3.190049551002597e-001, 9.085122161940182e-001);
  EXPECT_TRUE(Quat<double>::equal(Q2, Q2_expect));

  Quat<double> Q3(RotMat2);  // matrix --> quaternion
  Quat<double> Q3_expect(7.593187686640246e-002, 2.589898486876663e-001, 
    -3.190049551002597e-001, 9.085122161940182e-001);
  EXPECT_TRUE(Quat<double>::equal(Q3, Q3_expect));
}

TEST(Quaternion, Decomposition) {
  Vec3<double> axis(1.817168119886437e-001, 6.198030601001485e-001, 
    -7.634285604633715e-001);
  double angle = 8.621733203681206e-001;

  Mat4x4<double> RotMat1;
  // This is known correct due to a test previously done:
  Mat4x4<double>::rotateMatAxisAngle(RotMat1, axis, angle);
  Quat<double> Rot1(axis, angle);  // axis angle to quaternion

  Vec3<double> Scale(0.157613081677548f, 4.853756487228410f, 
    9.705927817606160f);
  Mat4x4<double> ScaleMat1;
  ScaleMat1.scaleMat(Scale[0], Scale[1], Scale[2]);

  Vec3<double> Trans(0.118997681558377f, 0.340385726666133f, 
    0.751267059305653f);
  Mat4x4<double> TransMat1;
  TransMat1.translationMat(Trans[0], Trans[1], Trans[2]);

  Mat4x4<double> SR, RS, TR, RT, ST, TS;
  Mat4x4<double>::mult(SR, ScaleMat1, RotMat1);
  Mat4x4<double>::mult(RS, RotMat1, ScaleMat1);
  Mat4x4<double>::mult(ST, ScaleMat1, TransMat1);
  Mat4x4<double>::mult(TS, TransMat1, ScaleMat1);
  Mat4x4<double>::mult(RT, RotMat1, TransMat1);
  Mat4x4<double>::mult(TR, TransMat1, RotMat1);

  Mat4x4<double> TSR, TRS, RST, RTS, STR, SRT;
  Mat4x4<double>::mult(TSR, TransMat1, SR);
  Mat4x4<double>::mult(TRS, TransMat1, RS);
  Mat4x4<double>::mult(RST, RotMat1, ST);
  Mat4x4<double>::mult(RTS, RotMat1, TS);
  Mat4x4<double>::mult(STR, ScaleMat1, TR);
  Mat4x4<double>::mult(SRT, ScaleMat1, RT);

  Vec3<double> Scale_ret;
  Quat<double> Rot_ret;
  Vec3<double> Trans_ret;

  Mat4x4<double> RS_ret, RST_ret, scale_mat_ret, trans_mat_ret;

  // NOTE: ROTATION IS ALWAYS CORRECT, scale should also be correct (thanks
  // to a hack of mine), but translation may not be correct.  HOWEVER, If
  // you do RST of the return matrix you will always get back the origional
  // matrix.

  Quat<double>::decompose(Trans_ret, Rot_ret, Scale_ret, TSR);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Quat<double>::equal(Rot_ret, Rot1));
  EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));

  Quat<double>::decompose(Trans_ret, Rot_ret, Scale_ret, TRS);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Quat<double>::equal(Rot_ret, Rot1));
  EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));

  Quat<double>::decompose(Trans_ret, Rot_ret, Scale_ret, RST);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Quat<double>::equal(Rot_ret, Rot1));
  // EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS

  Quat<double>::decompose(Trans_ret, Rot_ret, Scale_ret, RTS);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Quat<double>::equal(Rot_ret, Rot1));
  // EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS

  Quat<double>::decompose(Trans_ret, Rot_ret, Scale_ret, STR);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Quat<double>::equal(Rot_ret, Rot1));
  // EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS

  Quat<double>::decompose(Trans_ret, Rot_ret, Scale_ret, SRT);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Quat<double>::equal(Rot_ret, Rot1));
  // EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS
}
