# Nested Cages 

This code implements the following paper:
Leonardo Sacht and Etienne Vouga and Alec Jacobson: [Nested Cages] (http://www.cs.columbia.edu/cg/nested-cages/) ACM Transactions on Graphics, vol. 34, no. 6 (SIGGRAPH Asia 2015). 

It has been tested for Linux and Mac.

If you have any comments or questions, please contact
Leonardo Sacht by e-mail: leo@mtm.ufsc.br

------------- Getting Nested Cages -------------
git clone https://github.com/leokollersacht/nested_cages_demo.git

------------- Dependencies: -------------

1) libigl
git clone https://github.com/libigl/libigl.git
Note: This is a header library, it doesn't need to be compiled.

2) CGAL
Please install following the instructions at http://www.cgal.org/
Note: We are using CGAL 4.7 in this project, not guaranteed to be compatible with other versions.

3) Eigen
Download from http://eigen.tuxfamily.org/index.php?title=Main_Page
Note: this is a header-only library. For this project we are using Eigen 3.3.

4) Eltopo
git clone https://github.com/leokollersacht/eltopo.git
To compile, issue: 
cd eltopo3d/
make depend release
Note: This is our fork of Eltopo, contains small changes specific to our project

5) Tetgen
Please install following the instructions at http://wias-berlin.de/software/tetgen/

6) Meshfix
git clone https://github.com/evouga/collisiondetection.git
Notes: - This repository consists of the original Meshfix by Marco Attene 
with additional functionalities to support Eigen/libigl meshes. 
- You need to compile and link against OpenNL and SuperLU. For Linux systems
you should edit the Makefile to remove '-framework Accelerate' and set the 
paths to OpenNL and SuperLU accordingly.

7) Libcollisions (Implementation of Speculative Parallel Asynchronous Contact Mechanics - http://www.cs.columbia.edu/cg/spacm/spacm.html)
git clone https://github.com/evouga/collisiondetection.git
make

------------- Compiling Nested Cages: -------------

Once the above libraries are set (including changing their respecive cmake/findLIB file)
our program is compiled using the following commands:

cd nested_cages_demo/
cmake .
make

------------- Example usages -------------

./nested_cages_demo
Help information

./nested_cages_demo fertility.off 2 8000r 1000r None Volume test
Obtains 2 volume minimizing nested cages for fertility.off: one regular with 8000 faces and the other regular with 1000 faces. Output resulting cages to test_1.off and test_2.off

./nested_cages_demo fertility.off 2 8000 1000 None Volume test
The same as above, but outputs adaptive deimations (instead of regular)

./nested_cages_demo fertility.off 2 fertilitiy_8000.off fertilitiy_1000.off SurfARAP None test
Obtains 2 nested cages for fertility.off that minimize surface ARAP energy, using as input
decimations fertilitiy_8000.off and fertilitiy_1000.off