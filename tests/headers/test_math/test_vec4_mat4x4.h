//
//  test_vec2_mat2x2.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "test_unit/test_unit.h"
#include "icp/math/math_types.h"
#include "icp/math/math_base.h"

using icp::math::Vec4;
using icp::math::Mat4x4;

// Note: RotMat Axis angle is tested in test_quaternion

TEST(Vec4_Mat4x4, SimpleManipulation) {
  // Some practice variables
  Vec4<double> V_1(5.376671395461000e-001, 1.833885014595087e+000, 
    -2.258846861003648e+000, 8.621733203681206e-001);
  Vec4<double> V_2(3.187652398589808e-001, -1.307688296305273e+000, 
    -4.335920223056836e-001, 3.426244665386499e-001);
#ifdef ROW_MAJOR
  Mat4x4<double> M_1(3.578396939725761e+000, 7.254042249461056e-001, 
    -1.241443482163119e-001, 6.714971336080805e-001, 2.769437029884877e+000, 
    -6.305487318965619e-002, 1.489697607785465e+000, -1.207486922685038e+000,
    -1.349886940156521e+000, 7.147429038260958e-001, 1.409034489800479e+000, 
    7.172386513288385e-001, 3.034923466331855e+000, -2.049660582997746e-001, 
    1.417192413429614e+000, 1.630235289164729e+000);
  Mat4x4<double> M_2(4.888937703117894e-001, 2.938714670966581e-001, 
    -1.068870458168032e+000, 3.251905394561979e-001, 1.034693009917860e+000, 
    -7.872828037586376e-001, -8.094986944248755e-001, -7.549283191697034e-001,
    7.268851333832379e-001, 8.883956317576418e-001, -2.944284161994896e+000, 
    1.370298540095228e+000, -3.034409247860159e-001, -1.147070106969151e+000, 
    1.438380292815098e+000, -1.711516418853698e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat4x4<double> M_1(3.578396939725761e+000, 2.769437029884877e+000, 
    -1.349886940156521e+000, 3.034923466331855e+000, 7.254042249461056e-001,
    -6.305487318965619e-002, 7.147429038260958e-001, -2.049660582997746e-001,
    -1.241443482163119e-001, 1.489697607785465e+000, 1.409034489800479e+000,
    1.417192413429614e+000, 6.714971336080805e-001, -1.207486922685038e+000,
    7.172386513288385e-001, 1.630235289164729e+000);
  Mat4x4<double> M_2(4.888937703117894e-001, 1.034693009917860e+000, 
    7.268851333832379e-001, -3.034409247860159e-001, 2.938714670966581e-001, 
    -7.872828037586376e-001, 8.883956317576418e-001, -1.147070106969151e+000,
    -1.068870458168032e+000, -8.094986944248755e-001, -2.944284161994896e+000,
    1.438380292815098e+000, 3.251905394561979e-001, -7.549283191697034e-001,
    1.370298540095228e+000, -1.711516418853698e+000);
#endif

  double M_2_DET = Mat4x4<double>::det(M_2);
  EXPECT_APPROX_EQ(M_2_DET, -3.840132888329304e-001);

  M_2.transpose();
#ifdef ROW_MAJOR
  Mat4x4<double> M_2_expect(4.888937703117894e-001, 1.034693009917860e+000, 
    7.268851333832379e-001, -3.034409247860159e-001, 2.938714670966581e-001, 
    -7.872828037586376e-001, 8.883956317576418e-001, -1.147070106969151e+000,
    -1.068870458168032e+000, -8.094986944248755e-001, -2.944284161994896e+000, 
    1.438380292815098e+000, 3.251905394561979e-001, -7.549283191697034e-001, 
    1.370298540095228e+000, -1.711516418853698e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat4x4<double> M_2_expect(4.888937703117894e-001, 2.938714670966581e-001,
    -1.068870458168032e+000, 3.251905394561979e-001, 1.034693009917860e+000,
    -7.872828037586376e-001, -8.094986944248755e-001, -7.549283191697034e-001,
    7.268851333832379e-001, 8.883956317576418e-001, -2.944284161994896e+000, 
    1.370298540095228e+000, -3.034409247860159e-001, -1.147070106969151e+000,
    1.438380292815098e+000, -1.711516418853698e+000);
#endif
  EXPECT_TRUE(Mat4x4<double>::equal(M_2, M_2_expect));

  Mat4x4<double> M_3;
  Mat4x4<double>::mult(M_3, M_1, M_2);
#ifdef ROW_MAJOR
  Mat4x4<double> M_3_expect(2.313690316835965e+000, 2.725006513571085e+000, 
    4.531197261647137e+000, -3.245766731868780e+000, -6.495266052500589e-001, 
    2.620816957956779e+000, -4.083665709142190e+000, 3.441353201668658e+000,
    -1.722744886666602e+000, -3.641498727834494e+000, -3.512005130167420e+000, 
    3.889124363704155e-001, 4.388634886743703e-001, 9.236541541824801e-001, 
    8.524165717446142e-002, -1.437522370255467e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat4x4<double> M_3_expect(2.313690316835965e+000, -6.495266052500589e-001,
    -1.722744886666602e+000, 4.388634886743703e-001, 2.725006513571085e+000,
    2.620816957956779e+000, -3.641498727834494e+000, 9.236541541824801e-001,
    4.531197261647137e+000, -4.083665709142190e+000, -3.512005130167420e+000,
    8.524165717446142e-002, -3.245766731868780e+000, 3.441353201668658e+000,
    3.889124363704155e-001, -1.437522370255467e+000);
#endif
  EXPECT_TRUE(Mat4x4<double>::equal(M_3, M_3_expect));

  Vec4<double> V_3;
  Vec4<double>::add(V_3, V_1, V_2);
  Vec4<double> V_3_expect(8.564323794050808e-001, 5.261967182898131e-001, 
    -2.692438883309332e+000, 1.204797786906771e+000);
  EXPECT_TRUE(Vec4<double>::equal(V_3, V_3_expect));

  Vec4<double> V_4;
  Vec4<double>::sub(V_4, V_1, V_2);
  Vec4<double> V_4_expect(2.189018996871192e-001, 3.141573310900360e+000, 
    -1.825254838697965e+000, 5.195488538294706e-001);
  EXPECT_TRUE(Vec4<double>::equal(V_4, V_4_expect));

  Vec4<double> V_5;
  Vec4<double>::pairwiseMult(V_5, V_1, V_2);
  Vec4<double> V_5_expect(1.713895947017047e-001, -2.398149970355620e+000, 
    9.794179785414171e-001, 2.954016739549838e-001);
  EXPECT_TRUE(Vec4<double>::equal(V_5, V_5_expect));

  double VOT = Vec4<double>::dot(V_1, V_2);
  EXPECT_APPROX_EQ(VOT, -9.519407231575143e-001);

  Mat4x4<double> M_4;
  Mat4x4<double>::inverse(M_4, M_3);
#ifdef ROW_MAJOR
  Mat4x4<double> M_4_expect(4.224380293137507e+000, 1.756945738774290e+000,  
    3.299623525625032e+000, -4.439456166062089e+000, -1.072722341346062e+000, 
    -3.299182447285038e-001, -9.671415659374025e-001, 1.370627642625227e+000,
    -8.993283537885723e-001, -4.870271153233758e-001, -8.634340394938181e-001,
    6.310701746017395e-001, 5.470812859101863e-001, 2.955181703164883e-001, 
    3.347282993824935e-001, -1.132876865868168e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat4x4<double> M_4_expect(4.224380293137507e+000, -1.072722341346062e+000,
    -8.993283537885723e-001, 5.470812859101863e-001, 1.756945738774290e+000,
    -3.299182447285038e-001, -4.870271153233758e-001, 2.955181703164883e-001,
    3.299623525625032e+000, -9.671415659374025e-001, -8.634340394938181e-001,
    3.347282993824935e-001, -4.439456166062089e+000, 1.370627642625227e+000,
    6.310701746017395e-001, -1.132876865868168e+000);
#endif
  EXPECT_TRUE(Mat4x4<double>::equal(M_4, M_4_expect));
}

TEST(Vec4_Mat4x4, RotationMatrices) {
  // Correct value taken from quaternion test

  Vec3<double> axis(1.817168119886437e-001, 6.198030601001485e-001, 
    -7.634285604633715e-001);
  double angle = 8.621733203681206e-001;

  Mat4x4<double> RotMat1;
  Mat4x4<double>::rotateMatAxisAngle(RotMat1, axis, angle);

#ifdef ROW_MAJOR
  Mat4x4<double> RotMat1_expect(6.623201937964420e-001, 
    6.189709680704206e-001, 4.221455928650800e-001, 0, 
    -5.403086268696203e-001, 7.849403773940535e-001, -3.032081655673973e-001, 
    0, -5.190361727468607e-001, -2.726801464073089e-002, 
    8.543172167045705e-001, 0, 0, 0, 0, 1);
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat1, RotMat1_expect));
#endif
#ifdef COLUMN_MAJOR
  Mat4x4<double> RotMat1_expect(6.623201937964420e-001, 
    -5.403086268696203e-001, -5.190361727468607e-001, 0, 
    6.189709680704206e-001, 7.849403773940535e-001, -2.726801464073089e-002, 
    0, 4.221455928650800e-001,
    -3.032081655673973e-001, 8.543172167045705e-001, 0, 0, 0, 0, 1);
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat1, RotMat1_expect));
#endif

  // If the above is correct then use the general axis angle method to test
  // the specific axis tests
  axis.set(1, 0, 0);
  Mat4x4<double> RotMat2;
  Mat4x4<double>::rotateMatXAxis(RotMat2, angle);
  Mat4x4<double> RotMat2_expect;
  Mat4x4<double>::rotateMatAxisAngle(RotMat2_expect, axis, angle);
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat2, RotMat2_expect));

  axis.set(0, 1, 0);
  Mat4x4<double> RotMat3;
  Mat4x4<double>::rotateMatYAxis(RotMat3, -angle);
  Mat4x4<double> RotMat3_expect;
  Mat4x4<double>::rotateMatAxisAngle(RotMat3_expect, axis, -angle);
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat3, RotMat3_expect));

  axis.set(0, 0, 1);
  Mat4x4<double> RotMat4;
  Mat4x4<double>::rotateMatZAxis(RotMat4, 30*angle);
  Mat4x4<double> RotMat4_expect;
  Mat4x4<double>::rotateMatAxisAngle(RotMat4_expect, axis, 30*angle);
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat4, RotMat4_expect));

  // At least make sure the euler to mat --> mat to euler is invertable
  Vec3<double> euler(6.189709680704206e-001, 7.849403773940535e-001, 
    -2.726801464073089e-002);
  Mat4x4<double>::euler2RotMat(RotMat1, euler[0], euler[1], euler[2]);
  Vec3<double> euler_ret;
  Mat4x4<double>::rotMat2Euler(euler_ret[0], euler_ret[1], euler_ret[2],
    RotMat1);
  EXPECT_TRUE(Vec3<double>::equal(euler, euler_ret));
}

TEST(Vec4_Mat4x4, AffineInverses) {
  // The full 4x4 inverse is correct (if the above tests pass), use this as
  // the correct answer to test again

  Vec3<double> axis(1.817168119886437e-001, 6.198030601001485e-001, 
    -7.634285604633715e-001);
  double angle = 8.621733203681206e-001;

  Mat4x4<double> rot;
  Mat4x4<double>::rotateMatAxisAngle(rot, axis, angle);

  Mat4x4<double> trans;
  trans.translationMat(-9.671415659374025e-001, -2.726801464073089e-002,
    8.543172167045705e-001);

  Mat4x4<double> scale;
  scale.scaleMat(-3.032081655673973e-001, 2.955181703164883e-001,
    -4.439456166062089e+000);

  Mat4x4<double> RT;
  Mat4x4<double>::mult(RT, rot, trans);

  Mat4x4<double> RTS;
  Mat4x4<double>::mult(RTS, RT, scale);

  Mat4x4<double> AffineInverse;
  Mat4x4<double>::affineInverse(AffineInverse, RTS);

  Mat4x4<double> AffineInverse_expect;
  Mat4x4<double>::inverse(AffineInverse_expect, RTS);

  EXPECT_TRUE(Mat4x4<double>::equal(AffineInverse_expect, AffineInverse));

  Mat4x4<double> RT_inverse;
  Mat4x4<double>::affineRotationTranslationInverse(RT_inverse, RT);

  Mat4x4<double> RT_inverse_expect;
  Mat4x4<double>::inverse(RT_inverse_expect, RT);

  EXPECT_TRUE(Mat4x4<double>::equal(RT_inverse_expect, RT_inverse));

  Mat4x4<double> rot_inverse;
  Mat4x4<double>::affineRotationInverse(rot_inverse, rot);

  Mat4x4<double> rot_inverse_expect;
  Mat4x4<double>::inverse(rot_inverse_expect, rot);

  EXPECT_TRUE(Mat4x4<double>::equal(rot_inverse_expect, rot_inverse));
}
  
TEST(Vec4_Mat4x4, RightLeftMultTranslation) {
  // The full 4x4 inverse is correct (if the above tests pass), use this as
  // the correct answer to test again
    
  Vec3<double> axis(1.817168119886437e-001, 6.198030601001485e-001,
                    -7.634285604633715e-001);
  double angle = 8.621733203681206e-001;
    
  Mat4x4<double> rot;
  Mat4x4<double>::rotateMatAxisAngle(rot, axis, angle);
    
  Mat4x4<double> trans;
  Vec3<double> trans_vec(-9.671415659374025e-001, -2.726801464073089e-002,
                          8.543172167045705e-001);
  trans.translationMat(trans_vec[0], trans_vec[1], trans_vec[2]);
    
  Mat4x4<double> scale;
  scale.scaleMat(-3.032081655673973e-001, 2.955181703164883e-001,
                  -4.439456166062089e+000);
  {
    Mat4x4<double> RS;
    Mat4x4<double>::mult(RS, rot, scale);
      
    Mat4x4<double> TRS;
    Mat4x4<double>::mult(TRS, trans, RS);
      
    Mat4x4<double> RST;
    Mat4x4<double>::mult(RST, RS, trans);
      
    Mat4x4<double> TRS_quick;
    TRS_quick.set(RS);
    TRS_quick.leftMultTranslation(trans_vec[0], trans_vec[1], trans_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(TRS_quick, TRS));
      
    Mat4x4<double> RST_quick;
    RST_quick.set(RS);
    RST_quick.rightMultTranslation(trans_vec[0], trans_vec[1], trans_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(RST_quick, RST));
  }
    
  {
    Mat4x4<double> SR;
    Mat4x4<double>::mult(SR, scale, rot);
      
    Mat4x4<double> TSR;
    Mat4x4<double>::mult(TSR, trans, SR);
      
    Mat4x4<double> SRT;
    Mat4x4<double>::mult(SRT, SR, trans);
      
    Mat4x4<double> TSR_quick;
    TSR_quick.set(SR);
    TSR_quick.leftMultTranslation(trans_vec[0], trans_vec[1], trans_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(TSR_quick, TSR));
      
    Mat4x4<double> SRT_quick;
    SRT_quick.set(SR);
    SRT_quick.rightMultTranslation(trans_vec[0], trans_vec[1], trans_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(SRT_quick, SRT));
  }
}

TEST(Vec4_Mat4x4, RightLeftMultRotations) { 
#ifdef ROW_MAJOR
  Mat4x4<double> M_1(3.578396939725761e+000, 7.254042249461056e-001, 
    -1.241443482163119e-001, 6.714971336080805e-001, 2.769437029884877e+000, 
    -6.305487318965619e-002, 1.489697607785465e+000, -1.207486922685038e+000,
    -1.349886940156521e+000, 7.147429038260958e-001, 1.409034489800479e+000, 
    7.172386513288385e-001, 3.034923466331855e+000, -2.049660582997746e-001, 
    1.417192413429614e+000, 1.630235289164729e+000);
#endif
#ifdef COLUMN_MAJOR
  Mat4x4<double> M_1(3.578396939725761e+000, 2.769437029884877e+000, 
    -1.349886940156521e+000, 3.034923466331855e+000, 7.254042249461056e-001,
    -6.305487318965619e-002, 7.147429038260958e-001, -2.049660582997746e-001,
    -1.241443482163119e-001, 1.489697607785465e+000, 1.409034489800479e+000,
    1.417192413429614e+000, 6.714971336080805e-001, -1.207486922685038e+000,
    7.172386513288385e-001, 1.630235289164729e+000);
#endif
  double angle = 8.621733203681206e-001;
    
  Mat4x4<double> rot;
  Mat4x4<double> M_1_rot_inplace;
  Mat4x4<double> M_1_rot;

  // Left mult X-axis rotation
  Mat4x4<double>::rotateMatXAxis(rot, angle);
  M_1_rot_inplace.set(M_1);
  M_1_rot_inplace.leftMultRotateXAxis(angle);
  Mat4x4<double>::mult(M_1_rot, rot, M_1);
  EXPECT_TRUE(Mat4x4<double>::equal(M_1_rot, M_1_rot_inplace));

  // Left mult Y-axis rotation
  Mat4x4<double>::rotateMatYAxis(rot, angle);
  M_1_rot_inplace.set(M_1);
  M_1_rot_inplace.leftMultRotateYAxis(angle);
  Mat4x4<double>::mult(M_1_rot, rot, M_1);
  EXPECT_TRUE(Mat4x4<double>::equal(M_1_rot, M_1_rot_inplace));

  // Left mult Z-axis rotation
  Mat4x4<double>::rotateMatZAxis(rot, angle);
  M_1_rot_inplace.set(M_1);
  M_1_rot_inplace.leftMultRotateZAxis(angle);
  Mat4x4<double>::mult(M_1_rot, rot, M_1);
  EXPECT_TRUE(Mat4x4<double>::equal(M_1_rot, M_1_rot_inplace));

  // Right mult X-axis rotation
  Mat4x4<double>::rotateMatXAxis(rot, angle);
  M_1_rot_inplace.set(M_1);
  M_1_rot_inplace.rightMultRotateXAxis(angle);
  Mat4x4<double>::mult(M_1_rot, M_1, rot);
  EXPECT_TRUE(Mat4x4<double>::equal(M_1_rot, M_1_rot_inplace));

  // Right mult Y-axis rotation
  Mat4x4<double>::rotateMatYAxis(rot, angle);
  M_1_rot_inplace.set(M_1);
  M_1_rot_inplace.rightMultRotateYAxis(angle);
  Mat4x4<double>::mult(M_1_rot, M_1, rot);
  EXPECT_TRUE(Mat4x4<double>::equal(M_1_rot, M_1_rot_inplace));

  // Right mult Z-axis rotation
  Mat4x4<double>::rotateMatZAxis(rot, angle);
  M_1_rot_inplace.set(M_1);
  M_1_rot_inplace.rightMultRotateZAxis(angle);
  Mat4x4<double>::mult(M_1_rot, M_1, rot);
  EXPECT_TRUE(Mat4x4<double>::equal(M_1_rot, M_1_rot_inplace));
}
  
TEST(Vec4_Mat4x4, RightLeftMultScale) {
  // The full 4x4 inverse is correct (if the above tests pass), use this as
  // the correct answer to test again
    
  Vec3<double> axis(1.817168119886437e-001, 6.198030601001485e-001,
                    -7.634285604633715e-001);
  double angle = 8.621733203681206e-001;
    
  Mat4x4<double> rot;
  rot.rotateMatAxisAngle(axis, angle);
    
  Mat4x4<double> trans;
  trans.translationMat(-9.671415659374025e-001, -2.726801464073089e-002,
                        8.543172167045705e-001);
    
  Mat4x4<double> scale;
  Vec3<double> scale_vec(-3.032081655673973e-001, 2.955181703164883e-001,
                          -4.439456166062089e+000);
  scale.scaleMat(scale_vec[0], scale_vec[1], scale_vec[2]);
    
  {
    Mat4x4<double> RT;
    Mat4x4<double>::mult(RT, rot, trans);
      
    Mat4x4<double> RTS;
    Mat4x4<double>::mult(RTS, RT, scale);
      
    Mat4x4<double> SRT;
    Mat4x4<double>::mult(SRT, scale, RT);
      
    Mat4x4<double> RTS_quick;
    RTS_quick.set(RT);
    RTS_quick.rightMultScale(scale_vec[0], scale_vec[1], scale_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(RTS_quick, RTS));
      
    Mat4x4<double> SRT_quick;
    SRT_quick.set(RT);
    SRT_quick.leftMultScale(scale_vec[0], scale_vec[1], scale_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(SRT_quick, SRT));
  }
    
  {
    Mat4x4<double> TR;
    Mat4x4<double>::mult(TR, trans, rot);
      
    Mat4x4<double> TRS;
    Mat4x4<double>::mult(TRS, TR, scale);
      
    Mat4x4<double> STR;
    Mat4x4<double>::mult(STR, scale, TR);
      
    Mat4x4<double> TRS_quick;
    TRS_quick.set(TR);
    TRS_quick.rightMultScale(scale_vec[0], scale_vec[1], scale_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(TRS_quick, TRS));
      
    Mat4x4<double> STR_quick;
    STR_quick.set(TR);
    STR_quick.leftMultScale(scale_vec[0], scale_vec[1], scale_vec[2]);
    EXPECT_TRUE(Mat4x4<double>::equal(STR_quick, STR));
  }
}

TEST(Vec4_Mat4x4, Decomposition) {
  Vec3<double> axis(1.817168119886437e-001, 6.198030601001485e-001, 
    -7.634285604633715e-001);
  double angle = 8.621733203681206e-001;

  Mat4x4<double> RotMat1;
  // This is known correct due to a test previously done:
  RotMat1.rotateMatAxisAngle(axis, angle); 

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
  Mat4x4<double> RotMat_ret;
  Vec3<double> Trans_ret;

  Mat4x4<double> RS_ret, RST_ret, scale_mat_ret, trans_mat_ret;

  // NOTE: ROTATION IS ALWAYS CORRECT, scale should also be correct (thanks
  // to a hack of mine), but translation may not be correct.  HOWEVER, If
  // you do RST of the return matrix you will always get back the origional
  // matrix.

  Mat4x4<double>::decomposeRST(RotMat_ret, Scale_ret, Trans_ret, TSR);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat_ret, RotMat1));
  EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));

  Mat4x4<double>::decomposeRST(RotMat_ret, Scale_ret, Trans_ret, TRS);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat_ret, RotMat1));
  EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));

  Mat4x4<double>::decomposeRST(RotMat_ret, Scale_ret, Trans_ret, RST);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat_ret, RotMat1));
  // EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS

  Mat4x4<double>::decomposeRST(RotMat_ret, Scale_ret, Trans_ret, TS);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  // EXPECT_TRUE(Mat4x4<double>::equal(RotMat_ret, RotMat1));
  EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS

  Mat4x4<double>::decomposeRST(RotMat_ret, Scale_ret, Trans_ret, STR);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat_ret, RotMat1));
  // EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS

  Mat4x4<double>::decomposeRST(RotMat_ret, Scale_ret, Trans_ret, SRT);
  EXPECT_TRUE(Vec3<double>::equal(Scale_ret, Scale));
  EXPECT_TRUE(Mat4x4<double>::equal(RotMat_ret, RotMat1));
  // EXPECT_TRUE(Vec3<double>::equal(Trans_ret, Trans));  // Not true since T is not on LHS
}
