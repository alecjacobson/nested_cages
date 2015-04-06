#!/bin/bash
/*/../bin/ls > /dev/null
# BEGIN BASH SCRIPT
export PS4=""
set -o xtrace
TEMP="$0.cpp"
printf "//" | cat - $0 >$TEMP
#g++ -O3 -std=c++11 -fopenmp -o .main $TEMP -DNDEBUG -msse4.2 \
#clang++ -g -O0 -std=c++11 /opt/local/lib/gcc47/libstdc++.dylib -ferror-limit=4 -o .main $TEMP -msse4.2 \
clang++ -O3 -std=c++11 -o .main $TEMP -DNDEBUG -msse4.2 \
  -I. \
  -I/opt/local/include/eigen3/ -I/usr/local/igl/libigl/include \
  -I/usr/local/igl/libigl/external/AntTweakBar/include \
  -I/Applications/MATLAB_R2014b.app/extern/include/ \
  -I/usr/local/igl/libigl/external/AntTweakBar/src \
  -I/usr/local/igl/libigl/external/tetgen \
  -I/usr/local/igl/libigl/external/tinyxml2/ \
  -I/usr/local/igl/libigl/external/embree \
  -I/usr/local/igl/libigl/external/embree/embree \
  -L/usr/local/igl/libigl/external/embree/build -lembree -lsys \
  -L/usr/local/igl/libigl/external/tetgen -ltet \
  -I/usr/local/igl/libigl/external/glfw/include \
  -L/usr/local/igl/libigl/external/glfw/lib -lglfw3 -framework Carbon -framework QuartzCore -framework IOKit \
  -I/usr/local/igl/libigl/external/Singular_Value_Decomposition/ \
  -I/opt/local/include/ -I/usr/include/ \
  -I/usr/local/mosek/7/tools/platform/osx64x86/h \
  -L/usr/local/mosek/7/tools/platform/osx64x86/bin -lmosek64 \
  -L/opt/local/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt \
  -framework OpenGL \
  -framework GLUT \
  -framework AppKit \
  -L/opt/local/lib -lboost_thread-mt -lboost_system-mt \
  -L/usr/local/igl/libigl/external/AntTweakBar/lib -lAntTweakBar_clang \
  -L/Applications/MATLAB_R2014b.app/bin/maci64/ -lmx -lmat -lmex -leng \
  -L/opt/local/lib -lboost_program_options-mt -fno-math-errno && ./.main "$@"
#-DIGL_STATIC_LIBRARY -L/usr/local/igl/libigl/lib -liglviewer -ligl -liglmatlab -liglsvd3x3 -msse4.2 -fopenmp \
#rm -f .main
rm -f $TEMP
# END BASH SCRIPT
exit
*/

/* file exmaple.C

 This file is an example on how eiquadprog can be used
 by invoking solve_quadprog() function

 In order to compile this example, Eigen library must be installed
 on your system.
 
 The test problem is the following:
 
 Given:
 G =  2.1 0.0 1.0   g0^T = [6.0 1.0 1.0]
      1.5 2.2 0.0      
      1.2 1.3 3.1 
 Solve:
 min f(x) = 1/2 x G x + g0 x
 s.t.
   x_1 + 2*x_2 + x_3 = -4

   x_1 >= 0
   x_2 >= 0
   x_3 >= 0	
   -x_1 - x_2 >= -10
 
 The solution is x^T = [0 2 0] and f(x) = 6.4
 
 LICENSE
 
 Copyright (2010) Gael Guennebaud
 Copyright (2008) Angelo Furfaro
 


This file is a part of eiquadprog. 

uquadprog is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

uquadprog is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uquadprog; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/




#include <iostream>
#include <Eigen/Dense>

#include "eiquadprog.hpp"

using namespace Eigen;

template<typename Vec, typename Mat> void foo() {
  Mat G(3,3); 
  Vec g0(3);
  Mat CE(3,1);
  Vec ce0(1);
  Mat CI(3,4); 
  Vec ci0(4);
  Vec x(3);

  
  G(0,0)=2.1; G(0,1)=0.0; G(0,2)=1.0;
  G(1,0)=1.5; G(1,1)=2.2; G(1,2)=0.0;
  G(2,0)=1.2; G(2,1)=1.3; G(2,2)=3.1;
  
  
  g0(0)=6.0; g0(1)=1.0; g0(2)=1.0;

  CE(0,0)=1.0;  
  CE(1,0)=2.0;  
  CE(2,0)=-1.0; 
  
  ce0(0)=-4;

  CI(0,0)=1.0; CI(0,1)=0.0;CI(0,2)=0.0; CI(0,3)=-1.0;
  CI(1,0)=0.0; CI(1,1)=1.0;CI(1,2)=0.0; CI(1,3)=-1.0;
  CI(2,0)=0.0; CI(2,1)=0.0;CI(2,2)=1.0; CI(2,3)=0.0;


  ci0(0)=0.0; ci0(1)=0.0;ci0(2)=0.0; ci0(3)=10.0;


  std::cout << "f: " << solve_quadprog(G, g0,  CE, ce0,  CI, ci0, x) << std::endl;
  std::cout << "x: ";
  for (int i = 0; i < x.size(); i++)
    std::cout << x(i) << ' ';
  std::cout << std::endl;
}

int main(int argc, char** argv){
  	
  foo<Eigen::VectorXd, MatrixXd>();

}

