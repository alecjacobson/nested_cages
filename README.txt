# Nested Cages 

This code implements the following paper:
Leonardo Sacht and Etienne Vouga and Alec Jacobson: [Nested Cages] (http://www.cs.columbia.edu/cg/nested-cages/) ACM Transactions on Graphics, vol. 34, no. 6 (SIGGRAPH Asia 2015). 

It has been tested for Linux and Mac.

If you have any comments or questions, please contact
Leonardo Sacht by e-mail: leo@mtm.ufsc.br

----- Getting Nested Cages and submodules -----
git clone --recursive https://github.com/alecjacobson/nested_cages.git

------------- Dependencies: -------------

0) Libraries that are already included in this repository: Eltopo, Tetgen, Meshfix and Libcollisions
- Note: These libraries are submodules that are cloned and updated automatically (except for tetgen, that is included as a folder) if you clone our repository with the --recursive option. 
TO-DO: Cmake will compile and link them automatically.

1) libigl
git clone https://github.com/libigl/libigl.git
Notes: 
- This is a header library, it doesn't need to be compiled. 
- We recommend it to be located at /usr/local/include, so Cmake can find it automatically. You can also locate it at another path and edit cmake/FindLIBIGL accordingly. 

2) CGAL
Please install following the instructions at http://www.cgal.org/
Notes: 
- We are using CGAL 4.7 in this project, not guaranteed to be compatible with other versions.
- We recommend it to be installed at /usr/local/lib/cgal, so Cmake can find it automatically. You can also locate it at another path and edit cmake/FindCGAL accordingly. 


3) Eigen
Download from http://eigen.tuxfamily.org/index.php?title=Main_Page
Notes: 
- This is a header-only library. For this project we are using Eigen 3.3.
- We recommend it to be located at /usr/local/include/eigen3, so Cmake can find it automatically. You can also locate it at another path and edit cmake/FindEIGEN accordingly. 

------------- Compiling Nested Cages: -------------

Once the above libraries are set our program is compiled using the following commands:

cd nested_cages/
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