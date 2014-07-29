//
//  test_inverse_kinematics.cpp
//
//  Created by Jonathan Tompson on 7/15/14.
//
//  Test inverse kinematics code.  Here a simple 2-chain icp is simulated using
//  PSO as the optimizer.  Note that you can just as easily use BFGS or
//  PSOParallel (the API is more or less inter-changeable).  
//  Levenberg-Marquardt is problably not useful in this context.
//
//  For 
//

#include "test_unit/test_unit.h"
#include "test_unit/test_util.h"
#include "icp/math/icp.h"
#include "icp/clk/clk.h"
#include "icp/math/perlin_noise.h"
#include "icp/file_io/file_io.h"
#include "icp/math/common_optimization.h"

#include <sstream>
#include <thread>
#include <stdlib.h>
#include <GL/freeglut.h>

TEST(ICP, TEST_JACOBIAN) {
  double max_jacobian_err = 1e-3;
  EXPECT_LT(icp::math::ICP<double>::testJacobFunc(), max_jacobian_err);
}

// Some lazy globals
icp::math::ICP<float>* p_icp;
icp::clk::Clk clk;
icp::data_str::Vector<float> pc1, npc1;  // Point cloud 1 points
icp::data_str::Vector<float> pc2, npc2;  // Point cloud 2 points
double time1 = 0;
const int num_vertices = 1889;
const std::string bunny_vertices_filename = "bunny_vertices.bin";
const std::string bunny_normals_filename = "bunny_normals.bin";
int num_iterations = 29;
const int max_iterations = 30;
bool use_normals = false;
int win_width = 640;
int win_height = 640;
bool runtime_error_found = false;
int method = icp::math::ICPMethod::BFGS_ICP;

void RenderString(float x, float y, void *font, const char* str, 
  icp::math::Float3 const& rgb) {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, win_width, 0, win_height, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor3f(rgb[0], rgb[1], rgb[2]);
  glRasterPos2f(x,y);
  glutBitmapString(font, reinterpret_cast<const unsigned char *>(str));
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();   
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void renderScene(void) {
  float* pc2_transformed;
  if (num_iterations == 0) {
    pc2_transformed = &pc2[0];
  } else {
    icp::math::Float4x4 M_c;
    p_icp->num_iterations = num_iterations;
    p_icp->icp_method = (icp::math::ICPMethod)method;
    try {
      if (use_normals) {
        // Use Normals (better, but normals don't always exist)
        p_icp->match(M_c, &pc1[0], (pc1.size()/3), &pc2[0], (pc2.size()/3), M_c, 
          &npc1[0], &npc2[0]);
      } else {
        // Don't use normals (not as accurate... gets stuck in local minima)
        p_icp->match(M_c, &pc1[0], (pc1.size()/3), &pc2[0], (pc2.size()/3), M_c);
      }
    } catch (std::runtime_error& e) {
      std::cout << "exception caught!: " << e.what() << std::endl;
      runtime_error_found = true;
      glutLeaveMainLoop();
      return;
    }
    pc2_transformed = &p_icp->getPC2Transformed()[0];
  }
  float dist = 0;
  icp::math::Float3 delta;
  for (uint32_t i = 0; i < pc1.size()/3; i++) {
    delta.set(pc1[i*3+0] - pc2_transformed[i*3+0],
      pc1[i*3+1] - pc2_transformed[i*3+1],
      pc1[i*3+2] - pc2_transformed[i*3+2]);
    dist += delta.length();
  }

  // Render the points
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glLoadIdentity();
	// Set the camera
	gluLookAt(	0.0f, 0.0f,  2.0f,
			        0.0f, 0.0f,  0.0f,
			        0.0f, 1.0f,  0.0f );

  static float rot = 0.0f;
  double time2 = clk.getTime();
  rot += (float)(10 * (time2 - time1));
  time1 = time2;
  glRotatef(rot, 0, 1, 0);

  // Draw the 2 point clouds
  glPointSize(2.0f);
  glBegin(GL_POINTS); 
    glColor3f(1.0f, 1.0f, 1.0f);
    for (uint32_t i = 0; i < pc1.size()/3; i++) {
      glVertex3f(pc1[i*3+0], pc1[i*3+1], pc1[i*3+2]);
    }
    glColor3f(1.0f, 0.0f, 0.0f);
    for (uint32_t i = 0; i < pc2.size()/3; i++) {
      glVertex3f(pc2_transformed[i*3+0], pc2_transformed[i*3+1], 
        pc2_transformed[i*3+2]);
    }
  glEnd();

  // Draw the normals
  glBegin(GL_LINES); 
    glColor3f(0.0f, 0.0f, 1.0f);
    const float line_length = 0.02f;
    for (uint32_t i = 0; i < pc1.size()/3; i++) {
      glVertex3f(pc1[i*3+0], pc1[i*3+1], pc1[i*3+2]);
      glVertex3f(pc1[i*3+0] + line_length*npc1[i*3+0], 
                 pc1[i*3+1] + line_length*npc1[i*3+1], 
                 pc1[i*3+2] + line_length*npc1[i*3+2]);
    }
  glEnd();
  std::stringstream ss;
  ss << "num_iterations = " << num_iterations << ", use_normals(n) = ";
  ss << use_normals << ", icp_method(m) = " << method;
  RenderString(10, 10, GLUT_BITMAP_TIMES_ROMAN_24, ss.str().c_str(), 
    icp::math::Float3(1.0f, 1.0f, 1.0f));
  ss.str("");
  ss << "position err = " << dist;
  RenderString(10, 30, GLUT_BITMAP_TIMES_ROMAN_24, ss.str().c_str(), 
    icp::math::Float3(1.0f, 1.0f, 1.0f));

  glutSwapBuffers();
}

void changeSize(int w, int h) {
	if (h == 0) {
		h = 1;
  }
	float ratio =  w * 1.0f / h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(45.0f, ratio, 0.1f, 100.0f);
	glMatrixMode(GL_MODELVIEW);
  win_width = w;
  win_height = h;
}

void processSpecialKeys(int key, int x, int y) {
  switch(key) {
  case GLUT_KEY_RIGHT :
    num_iterations = (num_iterations + 1) % max_iterations;
    break;
  case GLUT_KEY_LEFT :
    num_iterations = num_iterations > 0 ? num_iterations-1 :
      (max_iterations-1);
    break;
  }
}

void processNormalKeys(unsigned char key, int x, int y) {
  switch (key) {
  case 27:  // ESCAPE
    glutLeaveMainLoop();
    break;
  case 'n':
    use_normals = !use_normals;
    break;
  case 'm':
    method = (method+1) % icp::math::ICPMethod::NUM_METHODS;
    break;
  }
}

TEST(ICP, SIMPLE_EXAMPLE) {
  p_icp = new icp::math::ICP<float>();
  p_icp->num_iterations = 1;
  p_icp->cos_normal_threshold = acosf((35.0f / 360.0f) * 2.0f * (float)M_PI);
  p_icp->min_distance_sq = 0.0001f;  // Avoid numerical issues
  p_icp->max_distance_sq = 100.0f;  // Helps avoid local minima
  p_icp->icp_method = icp::math::ICPMethod::BFGS_ICP;
  p_icp->verbose = false;

  // Load the Stanford bunny vertices (this will be our test point cloud)
  pc1.capacity(num_vertices * 3);  // Allocate enough space
  pc2.capacity(num_vertices * 3);
  npc1.capacity(num_vertices * 3);
  npc2.capacity(num_vertices * 3);
  pc1.resize(num_vertices * 3);  // Resize the vector manually
  pc2.resize(num_vertices * 3);
  npc1.resize(num_vertices * 3);
  npc2.resize(num_vertices * 3);
  icp::file_io::LoadArrayFromFile<float>(&pc1[0], num_vertices * 3, 
    bunny_vertices_filename);
  icp::file_io::LoadArrayFromFile<float>(&npc1[0], num_vertices * 3, 
    bunny_normals_filename);

  // Define a random rigid body transform:
  MERSINE_TWISTER_ENG eng;
  eng.seed(1000);
  UNIFORM_REAL_DISTRIBUTION dist_scale(0.8f, 1.2f);
  icp::math::Float3 scale(dist_scale(eng), dist_scale(eng), dist_scale(eng));
  UNIFORM_REAL_DISTRIBUTION dist_rot(-0.1f, +0.1f);
  icp::math::Float3 rot_euler(dist_rot(eng), dist_rot(eng), dist_rot(eng));
  UNIFORM_REAL_DISTRIBUTION dist_trans(-0.1f, +0.1f);
  icp::math::Float3 trans(dist_trans(eng), dist_trans(eng), dist_trans(eng));
  icp::math::Float4x4 M_pc2;

  icp::math::Float4x4::euler2RotMat(M_pc2, rot_euler[0], rot_euler[1],
    rot_euler[2]);
  M_pc2.leftMultScale(scale[0], scale[1], scale[2]);
  M_pc2.leftMultTranslation(trans[0], trans[1], trans[2]);

  std::cout << "Random transformation matrix for PC2:" << std::endl;
  M_pc2.printPrecise();

  // Create PC2 by transforming PC1 by the random rigid body transform
  icp::math::Float4x4 M_pc2_normal_mat;
  icp::math::Float4x4::inverse(M_pc2_normal_mat, M_pc2);
  M_pc2_normal_mat.transpose();
  icp::math::Float3 pt_pc1, pt_pc2, norm_pc1, norm_pc2;
  for (uint32_t i = 0; i < num_vertices; i++) {
    // Calculate the transformed position
    pt_pc1.set(pc1[i*3+0], pc1[i*3+1], pc1[i*3+2]);
    icp::math::Float3::affineTransformPos(pt_pc2, M_pc2, pt_pc1);
    pc2[i*3 + 0] = pt_pc2[0];
    pc2[i*3 + 1] = pt_pc2[1];
    pc2[i*3 + 2] = pt_pc2[2];
    
    // Calculate the transformed normal
    norm_pc1.set(npc1[i*3+0], npc1[i*3+1], npc1[i*3+2]);
    icp::math::Float3::affineTransformVec(norm_pc2, M_pc2_normal_mat, norm_pc1);
    norm_pc2.normalize();
    npc2[i*3 + 0] = norm_pc2[0];
    npc2[i*3 + 1] = norm_pc2[1];
    npc2[i*3 + 2] = norm_pc2[2];
  }

  // init GLUT and create Window
  int argc = 0;
  char ** argv = NULL;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("icp example");
	glutDisplayFunc(renderScene);
  glutReshapeFunc(changeSize);
  glutIdleFunc(renderScene);
  glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(processSpecialKeys);

	// Note: glutSetOption is only available with freeglut
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
                GLUT_ACTION_GLUTMAINLOOP_RETURNS);

  time1 = clk.getTime();
  glutMainLoop();

  EXPECT_FALSE(runtime_error_found);

  delete p_icp;
}