#!/bin/bash
/*/../bin/ls > /dev/null
# BEGIN BASH SCRIPT
export PS4=""
set -o xtrace
TEMP="$0.cpp"
printf "//" | cat - $0 >$TEMP
#g++ -O3 -std=c++11 -fopenmp -o .main $TEMP -DNDEBUG -msse4.2 \
#clang++ -O3 -std=c++11 -o .main $TEMP -DNDEBUG -msse4.2 \
clang++ -g -O0 -std=c++11 /opt/local/lib/gcc47/libstdc++.dylib -ferror-limit=4 -o .main $TEMP -msse4.2 \
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
#include <Eigen/Dense>
#include "linprog.h"
#include "mosek_linprog.h"
#include <iostream>
#include <stdexcept>


int main()
{
  using namespace Eigen;
  using namespace std;
  VectorXd c(3);
  c<< 0.0052293, -1.3553e-20, -0.008807;
  MatrixXd A(9,3);
  A<<0.000571443,-2.3183e-05,-0.000985265,0.000581547,-1.43161e-05,-0.000978401,0.000547457,-1.42998e-05,-0.000998138,0.000571443,2.3183e-05,-0.000985265,0.000581547,1.43161e-05,-0.000978401,0.000547457,1.42998e-05,-0.000998138,0.00060468,2.31543e-05,-0.000964086,0.000619,0,-0.000955257,0.00060468,-2.31543e-05,-0.000964086;
  VectorXd B(9);
  B<<0.00113896,0.00113802,0.00113824,0.00113896,0.00113802,0.00113824,0.001138,0.00113802,0.001138;
  {
    VectorXd x;
    linprog(c,-A,-B,MatrixXd(0,A.cols()),VectorXd(0,1),x);
    cout<<"x: "<<x.transpose()<<endl;
  }
  {
    VectorXd x;
    // B < A x
    mosek_linprog(c,A.sparseView(),B,VectorXd(),VectorXd(),VectorXd(),x);
    cout<<"x: "<<x.transpose()<<endl;
  }
}
