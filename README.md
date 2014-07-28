**icp - Jonathan's Iterative Closest Point Library**
---------
---------

**Overview**
--------

This library is a bunch of stuff ripped out of my private library (jtil) to implement the Iterative Closest Point algorithm.  Note: This library includes a variety of solutions for solving the transform for bringing two point clouds (PC1 and PC2 in to correspondence):

1. SVD to find the eigen vectors of the cross-covaraince matrix **(standard technique)**
2. A implementation of Shinji Umeyama's method described in "Least-Squares Estimation of Transformation Parameters Between Two Point Patterns" **(standard technique)**
3. A custom solution based on ICP and gradient descent **(non-standard technique)**

The third method *IS NOT* a generic ICP implementation, as I do not rely on any of the various closed form solutions for finding the 6DOF rigid body transform to bring two point clouds (PC1 and PC2) into correspondence.  Instead, I use BFGS to perform gradient descent so that the function to bring PC1 onto PC2 can be as non-linear as you want; i.e. you can parameterize shear, non-linear projection and any other second order effects.  

This formulation is particularly useful when dealing with real-world 3D scans (from the Kinect for instance), where you need to be able to compensate for non-linear depth and FOV mismatch between devices.  In general, I find that since the correspondance search (using KD trees) accounts for a huge percentage of the run-time, using BFGS is only slightly slower than the direct form quaternion methods, but is at the same time much more flexible.

For those who are curious, I use BFGS to solve the following objective function:

![Image of equation](eqn.png)

Where, `c` is the pose coefficient vector in `R^m`, `x_{PC1,i}` and `x_{PC2,i}` are the ith PC1 and PC2 correspondence pair in `R^3` and `M_c` is the user defined 4x3 linear transform for the given pose coefficient.

If you only want to solve for rotation and translation then use methods 1 and 2, however if you know that you need to solve for a global scale as well, then use method 3.  If you need shear and other non-linear components then you can modify the code of method 3.

**Compilation**
---------------

Building icp uses Visual Studio 2012 on Windows and there are no dependencies for the main library.  For the correspondence search I use the amazing nanoflann header library (included in the repo).   If you want to run the example you will need freeglut, however the pre-compiled x64 binaries are in this repo (so it the test project should just build and run).  The code is also cross platform (and compiles using gcc on Mac OS X), but I have not included any CMake files because I got sick of supporting them.

**Running**
---------------

If you build the test project (in icp.sln) then a simple glut window will show a the ICP algorithm progressing.  Press the right arrow to move forward to the next iteration.  For your own applications I would start from here.

**Style**
---------

This project follows the Google C++ style conventions: 

<http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml>
