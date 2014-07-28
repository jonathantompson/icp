package mathhelper;

import java.awt.*;
import java.awt.image.*;
import java.io.IOException;
import java.util.Hashtable;
import java.util.ArrayList;
import java.util.Arrays;

public class MathHelper {

	private final static double EPSILON = 0.00000001; // 1e-8
	
	// Note all matrices use OpenGL column-major form, ie:
	// M = |   M0    M4    M8    M12   |
	//     |   M1    M5    M9    M13   |
	//     |   M2    M6    M10   M14   |
	//     |   M3    M7    M11   M15   |
	
	public static void MatIdent(double[] matRet){
		matRet[0]  = 1.0;	matRet[4]  = 0.0;	matRet[8]  = 0.0;	matRet[12] = 0.0;
		matRet[1]  = 0.0;	matRet[5]  = 1.0;	matRet[9]  = 0.0;	matRet[13] = 0.0;
		matRet[2]  = 0.0;	matRet[6]  = 0.0;	matRet[10] = 1.0;	matRet[14] = 0.0;
		matRet[3]  = 0.0;	matRet[7]  = 0.0;	matRet[11] = 0.0;	matRet[15] = 1.0;
	}
	
	public static void CopyMat(double[] matDst, double[] matSrc){
		matDst[0]  = matSrc[0];		matDst[4]  = matSrc[4];		matDst[8]  = matSrc[8];		matDst[12]  = matSrc[12];
		matDst[1]  = matSrc[1];		matDst[5]  = matSrc[5];		matDst[9]  = matSrc[9];		matDst[13]  = matSrc[13];
		matDst[2]  = matSrc[2];		matDst[6]  = matSrc[6];		matDst[10]  = matSrc[10];	matDst[14]  = matSrc[14];
		matDst[3]  = matSrc[3];		matDst[7]  = matSrc[7];		matDst[11]  = matSrc[11];	matDst[15]  = matSrc[15];
	}
	
	public static void MatTranspose(double[] matDst, double[] matSrc){
		matDst[0]  = matSrc[0];		matDst[4]  = matSrc[1];		matDst[8]  = matSrc[2];		matDst[12]  = matSrc[3];
		matDst[1]  = matSrc[4];		matDst[5]  = matSrc[5];		matDst[9]  = matSrc[6];		matDst[13]  = matSrc[7];
		matDst[2]  = matSrc[8];		matDst[6]  = matSrc[9];		matDst[10]  = matSrc[10];	matDst[14]  = matSrc[11];
		matDst[3]  = matSrc[12];	matDst[7]  = matSrc[13];	matDst[11]  = matSrc[13];	matDst[15]  = matSrc[15];
	}
	
	public static void MatTranspose(double[] matRet){
		double temp;
		temp = matRet[1]; matRet[1] = matRet[4]; matRet[4] = temp;
		temp = matRet[2]; matRet[2] = matRet[8]; matRet[8] = temp;
		temp = matRet[3]; matRet[3] = matRet[12]; matRet[12] = temp;
		temp = matRet[6]; matRet[6] = matRet[9]; matRet[9] = temp;
		temp = matRet[7]; matRet[7] = matRet[13]; matRet[13] = temp;
		temp = matRet[11]; matRet[11] = matRet[14]; matRet[14] = temp;
	}
	
	public static void Mat3x3Transpose(double[] matRet){
		double temp;
		temp = matRet[1]; matRet[1] = matRet[3]; matRet[3] = temp;
		temp = matRet[2]; matRet[2] = matRet[6]; matRet[6] = temp;
		temp = matRet[5]; matRet[5] = matRet[7]; matRet[7] = temp;
	}
	
	public static void Mat3x3RotTranTo4x4(double[] matRet4x4, double[] matRot3x3, double[] vec3){
		matRet4x4[0]  = matRot3x3[0];	matRet4x4[4]  = matRot3x3[3];	matRet4x4[8]  = matRot3x3[6];	matRet4x4[12]  = vec3[0];
		matRet4x4[1]  = matRot3x3[1];	matRet4x4[5]  = matRot3x3[4];	matRet4x4[9]  = matRot3x3[7];	matRet4x4[13]  = vec3[1];
		matRet4x4[2]  = matRot3x3[2];	matRet4x4[6]  = matRot3x3[5];	matRet4x4[10] = matRot3x3[8];	matRet4x4[14]  = vec3[2];
		matRet4x4[3]  = 0.0;			matRet4x4[7]  = 0.0;			matRet4x4[11] = 0.0;			matRet4x4[15]  = 1.0;
	}
	
	public static void MatMult(double[] matRet, double[] matA, double[] matB){
		matRet[0]  = matA[0]*matB[0]   + matA[4]*matB[1]   + matA[8]*matB[2]   + matA[12]*matB[3];
		matRet[1]  = matA[1]*matB[0]   + matA[5]*matB[1]   + matA[9]*matB[2]   + matA[13]*matB[3];
		matRet[2]  = matA[2]*matB[0]   + matA[6]*matB[1]   + matA[10]*matB[2]  + matA[14]*matB[3];
		matRet[3]  = matA[3]*matB[0]   + matA[7]*matB[1]   + matA[11]*matB[2]  + matA[15]*matB[3];
		
		matRet[4]  = matA[0]*matB[4]   + matA[4]*matB[5]   + matA[8]*matB[6]   + matA[12]*matB[7];
		matRet[5]  = matA[1]*matB[4]   + matA[5]*matB[5]   + matA[9]*matB[6]   + matA[13]*matB[7];
		matRet[6]  = matA[2]*matB[4]   + matA[6]*matB[5]   + matA[10]*matB[6]  + matA[14]*matB[7];
		matRet[7]  = matA[3]*matB[4]   + matA[7]*matB[5]   + matA[11]*matB[6]  + matA[15]*matB[7];
		
		matRet[8]  = matA[0]*matB[8]   + matA[4]*matB[9]   + matA[8]*matB[10]  + matA[12]*matB[11];
		matRet[9]  = matA[1]*matB[8]   + matA[5]*matB[9]   + matA[9]*matB[10]  + matA[13]*matB[11];
		matRet[10] = matA[2]*matB[8]   + matA[6]*matB[9]   + matA[10]*matB[10] + matA[14]*matB[11];
		matRet[11] = matA[3]*matB[8]   + matA[7]*matB[9]   + matA[11]*matB[10] + matA[15]*matB[11];
		
		matRet[12] = matA[0]*matB[12]  + matA[4]*matB[13]  + matA[8]*matB[14]  + matA[12]*matB[15];
		matRet[13] = matA[1]*matB[12]  + matA[5]*matB[13]  + matA[9]*matB[14]  + matA[13]*matB[15];
		matRet[14] = matA[2]*matB[12]  + matA[6]*matB[13]  + matA[10]*matB[14] + matA[14]*matB[15];
		matRet[15] = matA[3]*matB[12]  + matA[7]*matB[13]  + matA[11]*matB[14] + matA[15]*matB[15];
	}
	
	// For 3D Coordinate: w = 1
	public static void Transform3DCoord(double[] vec4DRet, double[] matA, double[] vec3DB){
		vec4DRet[0]  = matA[0]*vec3DB[0]   + matA[4]*vec3DB[1]   + matA[8]*vec3DB[2]   + matA[12]; // vec3DB[3] = 1
		vec4DRet[1]  = matA[1]*vec3DB[0]   + matA[5]*vec3DB[1]   + matA[9]*vec3DB[2]   + matA[13];
		vec4DRet[2]  = matA[2]*vec3DB[0]   + matA[6]*vec3DB[1]   + matA[10]*vec3DB[2]  + matA[14];
		vec4DRet[3]  = matA[3]*vec3DB[0]   + matA[7]*vec3DB[1]   + matA[11]*vec3DB[2]  + matA[15];
	}
	public static void Transform3DCoordTo3D(double[] vec3DRet, double[] matA, double[] vec3DB){
		vec3DRet[0]  = matA[0]*vec3DB[0]   + matA[4]*vec3DB[1]   + matA[8]*vec3DB[2]   + matA[12]; // vec3DB[3] = 1
		vec3DRet[1]  = matA[1]*vec3DB[0]   + matA[5]*vec3DB[1]   + matA[9]*vec3DB[2]   + matA[13];
		vec3DRet[2]  = matA[2]*vec3DB[0]   + matA[6]*vec3DB[1]   + matA[10]*vec3DB[2]  + matA[14];
	}
	
	// For 3D Vector: w = 0
	public static void Transform3DVec(double[] vec4DRet, double[] matA, double[] vec3DB){
		vec4DRet[0]  = matA[0]*vec3DB[0]   + matA[4]*vec3DB[1]   + matA[8]*vec3DB[2]; // vec3DB[3] = 0
		vec4DRet[1]  = matA[1]*vec3DB[0]   + matA[5]*vec3DB[1]   + matA[9]*vec3DB[2];
		vec4DRet[2]  = matA[2]*vec3DB[0]   + matA[6]*vec3DB[1]   + matA[10]*vec3DB[2];
		vec4DRet[3]  = matA[3]*vec3DB[0]   + matA[7]*vec3DB[1]   + matA[11]*vec3DB[2];
	}
	
	public static void Transform4D(double[] vec3DRet, double[] matA, double[] vec3DB){
		vec3DRet[0]  = matA[0]*vec3DB[0]   + matA[4]*vec3DB[1]   + matA[8]*vec3DB[2]   + matA[12]*vec3DB[3];
		vec3DRet[1]  = matA[1]*vec3DB[0]   + matA[5]*vec3DB[1]   + matA[9]*vec3DB[2]   + matA[13]*vec3DB[3];
		vec3DRet[2]  = matA[2]*vec3DB[0]   + matA[6]*vec3DB[1]   + matA[10]*vec3DB[2]  + matA[14]*vec3DB[3];
		vec3DRet[3]  = matA[3]*vec3DB[0]   + matA[7]*vec3DB[1]   + matA[11]*vec3DB[2]  + matA[15]*vec3DB[3];
	}
	
	public static void CalcMatTranslation(double[] matRet, double[] vec3Trans) {
		matRet[0]  = 1.0;	matRet[4]  = 0.0;	matRet[8]  = 0.0;	matRet[12] = vec3Trans[0];
		matRet[1]  = 0.0;	matRet[5]  = 1.0;	matRet[9]  = 0.0;	matRet[13] = vec3Trans[1];
		matRet[2]  = 0.0;	matRet[6]  = 0.0;	matRet[10] = 1.0;	matRet[14] = vec3Trans[2];
		matRet[3]  = 0.0;	matRet[7]  = 0.0;	matRet[11] = 0.0;	matRet[15] = 1.0;
	}
	
	public static void CalcMatTranslation(double[] matRet, double dx, double dy, double dz) {
		matRet[0]  = 1.0;	matRet[4]  = 0.0;	matRet[8]  = 0.0;	matRet[12] = dx;
		matRet[1]  = 0.0;	matRet[5]  = 1.0;	matRet[9]  = 0.0;	matRet[13] = dy;
		matRet[2]  = 0.0;	matRet[6]  = 0.0;	matRet[10] = 1.0;	matRet[14] = dz;
		matRet[3]  = 0.0;	matRet[7]  = 0.0;	matRet[11] = 0.0;	matRet[15] = 1.0;
	}
		
	public static void CalcMatProj(double[] matRet, double fovy, double aspect, double zNear, double zFar ){
		// http://www.opengl.org/wiki/GluPerspective_code
		double ymax, xmax;
	    double temp, temp2, temp3, temp4;
	    ymax = zNear * Math.tan(fovy * Math.PI / 360.0);
	    //ymin = -ymax;
	    //xmin = -ymax * aspectRatio;
	    xmax = ymax * aspect;
	    CalcMatFrust(matRet, -xmax, xmax, -ymax, ymax, zNear, zFar);
	}
	
	public static void CalcMatFrust(double[] matRet, double left, double right, double bottom, double top, double znear, double zfar ){
		// http://www.opengl.org/wiki/GluPerspective_code
		double temp, temp2, temp3, temp4;
	    temp = 2.0 * znear;
	    temp2 = right - left;
	    temp3 = top - bottom;
	    temp4 = zfar - znear;
	    matRet[0] = temp / temp2;
	    matRet[1] = 0.0;
	    matRet[2] = 0.0;
	    matRet[3] = 0.0;
	    matRet[4] = 0.0;
	    matRet[5] = temp / temp3;
	    matRet[6] = 0.0;
	    matRet[7] = 0.0;
	    matRet[8] = (right + left) / temp2;
	    matRet[9] = (top + bottom) / temp3;
	    matRet[10] = (-zfar - znear) / temp4;
	    matRet[11] = -1.0;
	    matRet[12] = 0.0;
	    matRet[13] = 0.0;
	    matRet[14] = (-temp * zfar) / temp4;
	    matRet[15] = 0.0;
	}
	
	public static void CalcOrthogonalUpAndSide(double[] lookAtPt, double[] eyePt, double[] forward, double[] up, double[] side) {
        up[0] = 1.0; up[1] = 0.0; up[2] = 0.0; // default openGL up vector
		SubtractVec3D(forward,lookAtPt,eyePt);
		NormalizeVec3D(forward);
		CrossVec3D(side, forward, up);
		NormalizeVec3D(side);
		CrossVec3D(up, side, forward);
		NormalizeVec3D(up);
	}
	
	public static void CalcViewMat(double[] matView, double[] pos, double[] forward, double[] up, double[] side) {
		// http://www.cs.rutgers.edu/~decarlo/428/glu_man/lookat.html
		// http://www.opengl.org/wiki/GluLookAt_code

		// View matrix is inverse rotation 3x3 matrix (so transpose for rotation matrix), followed by translation
		double[] matRot = new double[16];
		matRot[0] = side[0];
		matRot[4] = side[1];
		matRot[8] = side[2];
		matRot[12] = 0.0;
		
		matRot[1] = up[0];
		matRot[5] = up[1];
		matRot[9] = up[2];
		matRot[13] = 0.0;
		
		matRot[2] = -forward[0];
		matRot[6] = -forward[1];
		matRot[10] = -forward[2];
		matRot[14] = 0.0;

		matRot[3] = 0.0; 
		matRot[7] = 0.0;
		matRot[11] = 0.0;
		matRot[15] = 1.0;
		
		double[] matTrans = new double[16];
		CalcMatTranslation(matTrans, -pos[0], -pos[1], -pos[2]);
		
		MatMult(matView, matRot, matTrans);
	}
	
	public static void SubtractVec3D(double[] vec3Ret, double[] vec3A, double[] vec3B) {
		vec3Ret[0] = vec3A[0] - vec3B[0];
		vec3Ret[1] = vec3A[1] - vec3B[1];
		vec3Ret[2] = vec3A[2] - vec3B[2];
	}
	
	public static void AddVec2D(double[] vec3Ret, double[] vec3A, double[] vec3B) {
		vec3Ret[0] = vec3A[0] + vec3B[0];
		vec3Ret[1] = vec3A[1] + vec3B[1];
	}
	
	public static void AddVec3D(double[] vec3Ret, double[] vec3A, double[] vec3B) {
		vec3Ret[0] = vec3A[0] + vec3B[0];
		vec3Ret[1] = vec3A[1] + vec3B[1];
		vec3Ret[2] = vec3A[2] + vec3B[2];
	}
	
	public static void AccumulateVec3D(double[] vec3Ret, double[] vec3A) {
		vec3Ret[0] += vec3A[0];
		vec3Ret[1] += vec3A[1];
		vec3Ret[2] += vec3A[2];
	}
	
	public static double DotVec3D(double[] vec3A, double[] vec3B) {
		return vec3A[0]*vec3B[0] + vec3A[1]*vec3B[1] + vec3A[2]*vec3B[2];
	}
	
	public static void CrossVec3D(double[] vec3Ret, double[] vec3A, double[] vec3B) {
		// From: http://en.wikipedia.org/wiki/Cross_product
		vec3Ret[0] = vec3A[1]*vec3B[2] - vec3A[2]*vec3B[1];
		vec3Ret[1] = vec3A[2]*vec3B[0] - vec3A[0]*vec3B[2];
		vec3Ret[2] = vec3A[0]*vec3B[1] - vec3A[1]*vec3B[0];
	}
	
	public static void SubtractVec2D(double[] vec2Ret, double[] vec2A, double[] vec2B) {
		vec2Ret[0] = vec2A[0] - vec2B[0];
		vec2Ret[1] = vec2A[1] - vec2B[1];
	}
	
	public static void SubtractVec2D(int[] vec2Ret, int[] vec2A, int[] vec2B) {
		vec2Ret[0] = vec2A[0] - vec2B[0];
		vec2Ret[1] = vec2A[1] - vec2B[1];
	}
	
	public static void CopyVec2D(double[] vec2DDest, double[] vec2DSrc) {
		vec2DDest[0] = vec2DSrc[0];
		vec2DDest[1] = vec2DSrc[1];
	}
	
	public static void CopyVec2D(int[] vec2DDest, int[] vec2DSrc) {
		vec2DDest[0] = vec2DSrc[0];
		vec2DDest[1] = vec2DSrc[1];
	}
	
	public static void CopyVec3D(double[] vec3DDest, double[] vec3DSrc) {
		vec3DDest[0] = vec3DSrc[0];
		vec3DDest[1] = vec3DSrc[1];
		vec3DDest[2] = vec3DSrc[2];
	}
	
	public static void CopyVec4D(double[] vec4DDest, double[] vec4DSrc) {
		vec4DDest[0] = vec4DSrc[0];
		vec4DDest[1] = vec4DSrc[1];
		vec4DDest[2] = vec4DSrc[2];
		vec4DDest[3] = vec4DSrc[3];
	}
	
	// http://www.cprogramming.com/tutorial/3d/rotationMatrices.html -> openGL is right-handed
	public static void CalcMatRotationXAxis(double[] matRet, double angle_rad) {
		  double cosA = Math.cos(angle_rad);
		  double sinA = Math.sin(angle_rad);
		  matRet[0] = 1.0;		matRet[4] = 0.0;	matRet[8] = 0.0;	matRet[12] = 0.0;
		  matRet[1] = 0.0;		matRet[5] = cosA;	matRet[9] = sinA;	matRet[13] = 0.0;
		  matRet[2] = 0.0;		matRet[6] = -sinA;	matRet[10] = cosA;	matRet[14] = 0.0;
		  matRet[3] = 0.0;		matRet[7] = 0.0;	matRet[0] = 0.0;	matRet[15] = 1.0;
	}
	
	public static void CalcMatRotationYAxis(double[] matRet, double angle_rad) {
		  double cosA = Math.cos(angle_rad);
		  double sinA = Math.sin(angle_rad);
		  matRet[0] = cosA;		matRet[4] = 0.0;	matRet[8] = -sinA;	matRet[12] = 0.0;
		  matRet[1] = 0.0;		matRet[5] = 1.0;	matRet[9] = 0.0;	matRet[13] = 0.0;
		  matRet[2] = sinA;		matRet[6] = 0.0;	matRet[10] = cosA;	matRet[14] = 0.0;
		  matRet[3] = 0.0;		matRet[7] = 0.0;	matRet[0] = 0.0;	matRet[15] = 1.0;
	}
	
	public static void CalcMatRotationZAxis(double[] matRet, double angle_rad) {
		  double cosA = Math.cos(angle_rad);
		  double sinA = Math.sin(angle_rad);
		  matRet[0] = cosA;		matRet[4] = sinA;	matRet[8] = 0.0;	matRet[12] = 0.0;
		  matRet[1] = -sinA;	matRet[5] = cosA;	matRet[9] = 0.0;	matRet[13] = 0.0;
		  matRet[2] = 0.0;		matRet[6] = 0.0;	matRet[10] =1.0;	matRet[14] = 0.0;
		  matRet[3] = 0.0;		matRet[7] = 0.0;	matRet[0] = 0.0;	matRet[15] = 1.0;
	}
	
	public static void CalcMatRotationAxisAngle(double[] matRet, double[] vec3Axis, double angle_rad) {
		// From": http://en.wikipedia.org/wiki/Rotation_matrix
		// pre-calculate to save computations
		double cosA = Math.cos(angle_rad);
		double sinA = Math.sin(angle_rad);
		double x = vec3Axis[0]; double y = vec3Axis[1]; double z = vec3Axis[2]; // A little slower, but easier to read
		double x2 = vec3Axis[0] * vec3Axis[0];
		double y2 = vec3Axis[1] * vec3Axis[1];
		double z2 = vec3Axis[2] * vec3Axis[2];
		matRet[0] = cosA + x2*(1-cosA);
		matRet[1] = y*x*(1-cosA) + z*sinA;
		matRet[2] = z*x*(1-cosA) - y*sinA;
		matRet[3] = 0.0;
		matRet[4] = x*y*(1-cosA) - z*sinA;
		matRet[5] = cosA + y2*(1-cosA);
		matRet[6] = z*y*(1-cosA) + x*sinA;
		matRet[7] = 0.0;
		matRet[8] = x*z*(1-cosA) + y*sinA;
		matRet[9] = y*z*(1-cosA) - x*sinA;
		matRet[10] = cosA + z2*(1-cosA);
		matRet[11] = 0.0;
		matRet[12] = 0.0;
		matRet[13] = 0.0;
		matRet[14] = 0.0;
		matRet[15] = 1.0;
	}
	
	public static void CalcMatRotationFromBasis(double[] matRet, double[] vec3A, double[] vec3B, double[] vec3C) {
		matRet[0]  = vec3A[0];	matRet[4]  = vec3B[0];	matRet[8]  = vec3C[0];	matRet[12] = 0.0;
		matRet[1]  = vec3A[1];	matRet[5]  = vec3B[1];	matRet[9]  = vec3C[1];	matRet[13] = 0.0;
		matRet[2]  = vec3A[2];	matRet[6]  = vec3B[2];	matRet[10] = vec3C[2];	matRet[14] = 0.0;
		matRet[3]  = 0.0;		matRet[7]  = 0.0;		matRet[11] = 0.0;		matRet[15] = 1.0;
	}
	
	public static void QuatToMat(double[] matRet, double[] quat){
		// Just to be save, normalize the quaternion
		NormalizeVec4D(quat);
		// Precalculate some common variables to save some time
		double x2 = quat[0] * quat[0]; 
		double y2 = quat[1] * quat[1]; 
		double z2 = quat[2] * quat[2]; 
		double xy = quat[0] * quat[1]; 
		double xz = quat[0] * quat[2];
		double yz = quat[1] * quat[2]; 
		double wx = quat[3] * quat[0]; 
		double wy = quat[3] * quat[1]; 
		double wz = quat[3] * quat[2];
		// From this point on, we MUST have a unit length quaternion
		matRet[0] = 1.0 - 2.0 * (y2 + z2); 
		matRet[1] = 2.0 * (xy - wz);
		matRet[2] = 2.0 * (xz + wy);
		matRet[3] = 0.0;
		matRet[4] = 2.0 * (xy + wz);
		matRet[5] = 1.0 - 2.0 * (x2 + z2);
		matRet[6] = 2.0 * (yz - wx);
		matRet[7] = 0.0;
		matRet[8] = 2.0 * (xz - wy); 
		matRet[9] = 2.0 * (yz + wx);
		matRet[10] = 1.0 - 2.0 * (x2 + y2);
		matRet[11] = 0.0; 
		matRet[12] = 0.0; 
		matRet[13] = 0.0; 
		matRet[14] = 0.0; 
		matRet[15] = 1.0;
	}
	
	public static void UnitQuaternion(double[] quatRet) {
		quatRet[0] = 0.0; // x
		quatRet[1] = 0.0; // y
		quatRet[2] = 0.0; // z
		quatRet[3] = 1.0; // w
	}
	
	public static void NormalizeVec4D(double[] vec) {
		double one_over_length = 1.0 / Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] + vec[3]*vec[3] + EPSILON);
		vec[0] *= one_over_length;
		vec[1] *= one_over_length;
		vec[2] *= one_over_length;
		vec[3] *= one_over_length;
	}
	
	public static void NormalizeVec3D(double[] vec) {
		double one_over_length = 1.0 / Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] + EPSILON);
		vec[0] *= one_over_length;
		vec[1] *= one_over_length;
		vec[2] *= one_over_length;
	}
	public static void NormalizeVec2D(double[] vec) {
		double one_over_length = 1.0 / Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + EPSILON);
		vec[0] *= one_over_length;
		vec[1] *= one_over_length;;
	}
	
	public static double LengthVec3D(double[] vec) {
		return Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	}
	
	public static void MultVec3DScalar(double[] vec, double scalar) {
		vec[0] *= scalar;
		vec[1] *= scalar;
		vec[2] *= scalar;
	}
	
	public static void MultVec2DScalar(double[] vec, double scalar) {
		vec[0] *= scalar;
		vec[1] *= scalar;
	}
	
	public static void PrintVec3D(double[] vec) {
		System.out.println("<" + vec[0] + ", " + vec[1] + ", " + vec[2] + ">");
	}
	
	public static void SaveMatToFile(double[] mat, String filename) {
		try {
			java.io.FileOutputStream fos = new java.io.FileOutputStream(filename);
			java.io.ObjectOutputStream outStream = new java.io.ObjectOutputStream(fos);
			for(int i = 0; i < mat.length; i ++)
				outStream.writeDouble(mat[i]);
			outStream.flush(); outStream.close(); fos.close();
		} catch(IOException e) {
			System.out.println("SaveMatToFile: Caught IOException when trying to write to file: " + e.getMessage());
			System.exit(-1);
		}
	}
	
	public static void LoadMatFromFile(double[] mat, String filename) {
		try {
			java.io.FileInputStream fis = new java.io.FileInputStream(filename);
			java.io.ObjectInputStream inStream = new java.io.ObjectInputStream(fis);
			for(int i = 0; i < mat.length; i ++)
				mat[i] = inStream.readDouble();
			inStream.close(); inStream.close();
		} catch(IOException e) {
			System.out.println("SaveMatToFile: Caught IOException when trying to write to file: " + e.getMessage());
			System.exit(-1);
		}
	}
	
	public static void JoinCoordSystems(double[][][] corresp, double[] mat1, double[] mat2) { // corresp[curCalKinect][curCalCor][0] 
		double[] pA_k1 = new double[3];
		double[] pB_k1 = new double[3];
		double[] pC_k1 = new double[3];
		double[] pD_k1 = new double[3];
		double[] pA_k2 = new double[3];
		double[] pB_k2 = new double[3];
		double[] pC_k2 = new double[3];
		double[] pD_k2 = new double[3];
		
		// TRANSFORM MAT2's inputs, but leave mat1 transformed by identity matrix
		CopyVec3D(pA_k1, corresp[0][0]);
		CopyVec3D(pB_k1, corresp[0][1]);
		CopyVec3D(pC_k1, corresp[0][2]);
		CopyVec3D(pD_k1, corresp[0][3]);
		Transform3DCoordTo3D(pA_k2, mat2, corresp[1][0]);
		Transform3DCoordTo3D(pB_k2, mat2, corresp[1][1]);
		Transform3DCoordTo3D(pC_k2, mat2, corresp[1][2]);
		Transform3DCoordTo3D(pD_k2, mat2, corresp[1][3]);
		
		// Now calculate an orthonormal basis for each kinect using these 4 points (only 3 points is actually needed)
		double[] v1 = new double[3];
		double[] v2 = new double[3];
		double[] v3 = new double[3];
		// Find kinect 1's rotation matrix
		SubtractVec3D(v1, pA_k1, pB_k1 ); NormalizeVec3D(v1);
		SubtractVec3D(v2, pC_k1, pB_k1 ); NormalizeVec3D(v2);
		CrossVec3D(v3, v1, v2); NormalizeVec3D(v3); // Orthogonal to V1, V2
		CrossVec3D(v2, v3, v1); NormalizeVec3D(v3); // Now cross again to make basis
		double[] RotA = new double[16];
		CalcMatRotationFromBasis(RotA, v1, v2, v3);
		// Find kinect 2's rotation matrix
		SubtractVec3D(v1, pA_k2, pB_k2 ); NormalizeVec3D(v1);
		SubtractVec3D(v2, pC_k2, pB_k2 ); NormalizeVec3D(v2);
		CrossVec3D(v3, v1, v2); NormalizeVec3D(v3); // Orthogonal to V1, V2
		CrossVec3D(v2, v3, v1); NormalizeVec3D(v3); // Now cross again to make basis
		double[] RotB = new double[16];
		CalcMatRotationFromBasis(RotB, v1, v2, v3);
		
		// Now bring A into B's coordinate frame: R = B * A^-1;
		MatTranspose(RotA); // Rotation matrix is orthonormal, so just take transpose to find inverse
		MatMult(mat1, RotB, RotA ); // RotA^-1 first, followed by RotB
		
		// Now get the translation by bringing one of the points into A's coordinate system
		// To help with accuracy, do this for all 4 points and take the average
		Transform3DCoordTo3D(pA_k1, mat1, corresp[0][0]); SubtractVec3D(v1, pA_k2, pA_k1 );
		CopyVec3D(v2, v1);
		Transform3DCoordTo3D(pB_k1, mat1, corresp[0][1]); SubtractVec3D(v1, pB_k2, pB_k1 );
		AccumulateVec3D(v2, v1);
		Transform3DCoordTo3D(pC_k1, mat1, corresp[0][2]); SubtractVec3D(v1, pC_k2, pC_k1 );
		AccumulateVec3D(v2, v1);
		Transform3DCoordTo3D(pD_k1, mat1, corresp[0][3]); SubtractVec3D(v1, pD_k2, pD_k1 );
		AccumulateVec3D(v2, v1);
		MultVec3DScalar(v2, 1.0/4.0);
		CalcMatTranslation(RotA, v2); // use RotA as a temporary matrix
		
		// Now calculate composite transformation
		MatMult(RotB, RotA, mat1); // mat1 (rotation) first, followed by RotA (translation)54
		CopyMat(mat1, RotB);
		
		// Let's check some values:
		System.out.println("JoinCoordSystems results: ");
		Transform3DCoordTo3D(pA_k1, mat1, corresp[0][0]);
		Transform3DCoordTo3D(pA_k2, mat2, corresp[1][0]);
		System.out.print("mat1 * corresp[0][0] = "); PrintVec3D(pA_k1);
		System.out.print("mat2 * corresp[1][0] = "); PrintVec3D(pA_k2);
		Transform3DCoordTo3D(pB_k1, mat1, corresp[0][1]);
		Transform3DCoordTo3D(pB_k2, mat2, corresp[1][1]);
		System.out.print("mat1 * corresp[0][1] = "); PrintVec3D(pB_k1);
		System.out.print("mat2 * corresp[1][1] = "); PrintVec3D(pB_k2);
		Transform3DCoordTo3D(pC_k1, mat1, corresp[0][2]);
		Transform3DCoordTo3D(pC_k2, mat2, corresp[1][2]);
		System.out.print("mat1 * corresp[0][2] = "); PrintVec3D(pC_k1);
		System.out.print("mat2 * corresp[1][2] = "); PrintVec3D(pC_k2);
		Transform3DCoordTo3D(pD_k1, mat1, corresp[0][3]);
		Transform3DCoordTo3D(pD_k2, mat2, corresp[1][3]);
		System.out.print("mat1 * corresp[0][3] = "); PrintVec3D(pD_k1);
		System.out.print("mat2 * corresp[1][3] = "); PrintVec3D(pD_k2);
		
	}
	
	// NEED 4 POINTS SINCE THERE ARE 12 LINEAR SYSTEMS OF EQUATIONS (4X4 MATRIX)
	// THIS IS A VERY GENERAL FUNCTION THAT JOINS TWO ARBITRARY COORDINATE SYSTEMS
	public static void JoinCoordSystemsGeneral(double[][][] corresp, double[] mat1, double[] mat2) { // corresp[curCalKinect][curCalCor][0] 
		double[] pA_k1 = new double[3];
		double[] pB_k1 = new double[3];
		double[] pC_k1 = new double[3];
		double[] pD_k1 = new double[3];
		double[] pA_k2 = new double[3];
		double[] pB_k2 = new double[3];
		double[] pC_k2 = new double[3];
		double[] pD_k2 = new double[3];
		
		// TRANSFORM MAT2's inputs, but leave mat1 transformed by identity matrix
		CopyVec3D(pA_k1, corresp[0][0]);
		CopyVec3D(pB_k1, corresp[0][1]);
		CopyVec3D(pC_k1, corresp[0][2]);
		CopyVec3D(pD_k1, corresp[0][3]);
		Transform3DCoordTo3D(pA_k2, mat2, corresp[1][0]);
		Transform3DCoordTo3D(pB_k2, mat2, corresp[1][1]);
		Transform3DCoordTo3D(pC_k2, mat2, corresp[1][2]);
		Transform3DCoordTo3D(pD_k2, mat2, corresp[1][3]);
				
//		For a system of 4 points in 2 GENERAL coordinate systems (they don't even need to be orthogonal)
//		point A:   (ax1,ay1,az1)   and   (ax2,ay2,az2)	  
//		point B:   (bx1,by1,bz1)   and   (bx2,by2,bz2)
//		point C:   (cx1,cy1,cz1)   and   (cx2,cy2,cz2)
//		point D:   (dx1,dy1,dz1)   and   (dx2,dy2,dz2)
//		Set up:  (ax2,ay2,az2) = M * (ax1,ay1,az1) + O, Where M is 3x3 and O is translation
//		         (bx2,by2,bz2) = M * (bx1,by1,bz1) + O
//               (cx2,cy2,cz2) = M * (cx1,cy1,cz1) + O
//               (dx2,dy2,dz2) = M * (dx1,dy1,dz1) + O
//      --> 12 EQUATIONS (each above is split into 3 scalars) + 12 UNKNOWNS (3X3 + 3)
//      --> Results in a linear system of P * M = Q --> M = P^-1 * Q
		
//	       | (bx1-ax1) (by1-ay1) (bz1-az1) |
//	   P = | (cx1-ax1) (cy1-ay1) (cz1-az1) |
//	       | (dx1-ax1) (dy1-ay1) (dz1-az1) |
		double[] P = new double[16];
		MatIdent(P);
		P[0] = (pB_k1[0] - pA_k1[0]);	P[4] = (pB_k1[1] - pA_k1[1]);	P[8] = (pB_k1[2] - pA_k1[2]);
		P[1] = (pC_k1[0] - pA_k1[0]);	P[5] = (pC_k1[1] - pA_k1[1]);	P[9] = (pC_k1[2] - pA_k1[2]);
		P[2] = (pD_k1[0] - pA_k1[0]);	P[6] = (pD_k1[1] - pA_k1[1]);	P[10] = (pD_k1[2] - pA_k1[2]);

//	       | (bx2-ax2) (by2-ay2) (bz2-az2) |
//	   Q = | (cx2-ax2) (cy2-ay2) (cz2-az2) |
//	       | (dx2-ax2) (dy2-ay2) (dz2-az2) |
		double[] Q = new double[16];
		MatIdent(Q);
		Q[0] = (pB_k2[0] - pA_k2[0]);	Q[4] = (pB_k2[1] - pA_k2[1]);	Q[8] = (pB_k2[2] - pA_k2[2]);
		Q[1] = (pC_k2[0] - pA_k2[0]);	Q[5] = (pC_k2[1] - pA_k2[1]);	Q[9] = (pC_k2[2] - pA_k2[2]);
		Q[2] = (pD_k2[0] - pA_k2[0]);	Q[6] = (pD_k2[1] - pA_k2[1]);	Q[10] = (pD_k2[2] - pA_k2[2]);
		
//	      P * M  = Q
//        M  = P^-1 * Q
		double[] Pinv = new double[16];
		gluInvertMatrix(Pinv, P);
		MatMult(P, Pinv, Q); // USE P AS A TEMPORARY MATRIX
		MatTranspose(P);
		
//		(Ox,Oy,Oz)  =  (ax2,ay2,az2)  -  (ax1,ay1,az1) * M
		double[] vecTemp = new double[3];
		double[] trans = new double[3];
		Transform3DCoordTo3D(vecTemp, P, pB_k1);
		SubtractVec3D(trans, pB_k2, vecTemp);
		
		P[12] = trans[0];
		P[13] = trans[1];
		P[14] = trans[2];
		
		CopyMat(mat1, P); // USE PINV AS A TEMPORARY MATRIX
		
		// Let's check some values:
		System.out.println("JoinCoordSystems results: ");
		Transform3DCoordTo3D(pA_k1, mat1, corresp[0][0]);
		Transform3DCoordTo3D(pA_k2, mat2, corresp[1][0]);
		System.out.print("mat1 * corresp[0][0] = "); PrintVec3D(pA_k1);
		System.out.print("mat2 * corresp[1][0] = "); PrintVec3D(pA_k2);
		Transform3DCoordTo3D(pB_k1, mat1, corresp[0][1]);
		Transform3DCoordTo3D(pB_k2, mat2, corresp[1][1]);
		System.out.print("mat1 * corresp[0][1] = "); PrintVec3D(pB_k1);
		System.out.print("mat2 * corresp[1][1] = "); PrintVec3D(pB_k2);
		Transform3DCoordTo3D(pC_k1, mat1, corresp[0][2]);
		Transform3DCoordTo3D(pC_k2, mat2, corresp[1][2]);
		System.out.print("mat1 * corresp[0][2] = "); PrintVec3D(pC_k1);
		System.out.print("mat2 * corresp[1][2] = "); PrintVec3D(pC_k2);
		Transform3DCoordTo3D(pD_k1, mat1, corresp[0][3]);
		Transform3DCoordTo3D(pD_k2, mat2, corresp[1][3]);
		System.out.print("mat1 * corresp[0][3] = "); PrintVec3D(pD_k1);
		System.out.print("mat2 * corresp[1][3] = "); PrintVec3D(pD_k2);
		
	}
	
	// 4x4 matrix invert stolen from here: http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
	// This is a direct form inverse and avoid gauss-jordan
	public static Boolean gluInvertMatrix(double[] inv, double[] m) {
		
		inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15] + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
		inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15] - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
		inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15] + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
		inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14] - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
		inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15] - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
		inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15] + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
		inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15] - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
		inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14] + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
		inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15] + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
		inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15] - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
		inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15] + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
		inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14] - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
		inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11] - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
		inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11] + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
		inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11] - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
		inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10] + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];
	
		double det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
		if (Math.abs(det) < 0.00001) {
			System.out.println("Warning: determinate in gluInvertMatrix is < 0.00001: " + det);
			return false;
		}
	
		det = 1.0 / det;
	
		for (int i = 0; i < 16; i++)
		        inv[i] *= det;
	
		return true;
	}
	
	public static Boolean MatEqual(double[] matA, double[] matB, double epsilon) {
		for (int i = 0; i < matA.length; i ++ )
			if(Math.abs(matA[i] - matB[i]) > epsilon)
				return false;
		return true;
	}
	
	public static int maxArray( int[] arr ) {
		int retVal = arr[0];
		for (int i = 1; i < arr.length; i ++) {
			if(arr[i] > retVal)
				retVal = arr[i];
		}
		return retVal;
	}
	
	public static int minArray( int[] arr ) {
		int retVal = arr[0];
		for (int i = 1; i < arr.length; i ++) {
			if(arr[i] < retVal)
				retVal = arr[i];
		}
		return retVal;
	}
	
	// Quick test structures to compare against Matlab output (or here http://bmanolov.free.fr/matrixcalc.php)
	public static void TestMath() {
		double[] matA = new double[16];
		double[] matB = new double[16];
		double[] matC = new double[16];
		double[] matD = new double[16];
		double[] vecA = new double[4];
		double[] vecB = new double[4];
		double[] vecC = new double[4];
		double doubleA;
		
		matA[0] = 1;   matA[4] = 2;   matA[8] = 3;   matA[12] = 4;
		matA[1] = 5;   matA[5] = 6;   matA[9] = 7;   matA[13] = 8;
		matA[2] = 9;   matA[6] = 10;  matA[10] = 11; matA[14] = 12;
		matA[3] = 13;  matA[7] = 14;  matA[11] = 15; matA[15] = 16;
		
		matB[0] = 17;  matB[4] = 18;  matB[8] = 19;  matB[12] = 20;
		matB[1] = 21;  matB[5] = 22;  matB[9] = 23;  matB[13] = 24;
		matB[2] = 25;  matB[6] = 26;  matB[10] = 27; matB[14] = 28;
		matB[3] = 29;  matB[7] = 30;  matB[11] = 31; matB[15] = 32;
		
		vecA[0] = 33;  vecA[1] = 34;  vecA[2] = 35;  vecA[3] = 36;
		
		vecB[0] = 37;  vecB[1] = 38;  vecB[2] = 39;  vecB[3] = 40;
		
		MatMult(matC, matA, matB);
		// expect C =
		// 250    260    270    280   
		// 618    644    670    696   
		// 986   1028   1070   1112   
		//1354   1412   1470   1528   
		
		Transform4D(vecC, matA, vecA);
		// expect C =
		//  350   902   1454   2006  
		
		Transform3DVec(vecC, matA, vecA);
		// expect C =
		//  206   614   1022   1430    
		
		Transform3DCoord(vecC, matA, vecA);
		// expect C =
		//  210   622   1034   1446   
		
		// http://members.tripod.com/c_carleton/dotprod.html
		doubleA = DotVec3D(vecA, vecB);
		// expect doubleA = 3878
		
		// http://www.analyzemath.com/vector_calculators/vector_cross_product.html
		CrossVec3D(vecC, vecA, vecB);
		// expect vecC =
		//   -4     8     -4
		
		double[][][] corresp = new double[2][4][3];
		//        1st "known" frame        
		corresp[0][0][0] =  0.1; corresp[0][0][1] =  1.2; corresp[0][0][2] =  3.5;
		corresp[0][1][0] =  2.1; corresp[0][1][1] =  1.1; corresp[0][1][2] =  2.0;
		corresp[0][2][0] =  1.8; corresp[0][2][1] =  0.2; corresp[0][2][2] =  0.9;
		corresp[0][3][0] =  3.1; corresp[0][3][1] =  1.1; corresp[0][3][2] =  0.0;
		//        2nd "unknown" frame
		corresp[1][0][0] =  1.2; corresp[1][0][1] =  0.6; corresp[1][0][2] = -0.8;
		corresp[1][1][0] =  1.0; corresp[1][1][1] =  1.0; corresp[1][1][2] = -0.5;
		corresp[1][2][0] = -0.3; corresp[1][2][1] =  1.5; corresp[1][2][2] =  1.3;
		corresp[1][3][0] = -1.0; corresp[1][3][1] =  0.2; corresp[1][3][2] =  2.3;
		MatIdent(matA); MatIdent(matB); 
		// JoinCoordSystems(corresp, matA, matB);
		
		// Another way to test it: just make up a matrix and see if the function can give us that matrix back.
		vecA[0] = 1.0; vecA[1] = 1.0; vecA[1] = 1.0; NormalizeVec3D(vecA);
		CalcMatRotationAxisAngle(matC, vecA, 2.0 * Math.PI * 31.7 / 360.0); // Some rotation
		matC[12] = 1.2; matC[13] = -9.6; matC[14] = -2.5; // Some translation
		Transform3DCoordTo3D(vecA, matC, corresp[1][0]); CopyVec3D(corresp[0][0], vecA); // pA_k1 = matC * pA_k2
		Transform3DCoordTo3D(vecA, matC, corresp[1][1]); CopyVec3D(corresp[0][1], vecA); // pB_k1 = matC * pB_k2
		Transform3DCoordTo3D(vecA, matC, corresp[1][2]); CopyVec3D(corresp[0][2], vecA); // pC_k1 = matC * pC_k2
		Transform3DCoordTo3D(vecA, matC, corresp[1][3]); CopyVec3D(corresp[0][3], vecA); // pD_k1 = matC * pD_k2
		
		gluInvertMatrix(matD, matC);
		JoinCoordSystems(corresp, matA, matB); // Now matA should equal matC!!
		if(MatEqual(matA,matD, 0.000001))
			System.out.println("MatA == MatB is true");
		
	}
	
}
