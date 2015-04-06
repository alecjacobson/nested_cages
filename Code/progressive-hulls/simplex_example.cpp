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
#include "SimplexSolver.h"
#include "SimplexSolver.cpp"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>


int main()
{
  using namespace Eigen;
  using namespace std;
	SimplexSolver *solver1 = NULL, *solver2 = NULL;
	MatrixXd constraints(3, 3);
	VectorXd objectiveFunction(2);
	
	try {
		/*
			Maximization problem
		*/
		objectiveFunction <<	1,
								2;

		constraints <<		2,	3,	34,
							1,	5,	45,
							1,	0,	15;

		
		solver1 = new SimplexSolver(SIMPLEX_MAXIMIZE, objectiveFunction, constraints);

		if (solver1->hasSolution()) {
			cout << "The maximum is: " << solver1->getOptimum() << endl;
			cout << "The solution is: " << solver1->getSolution().transpose() << endl;
		} else {
			cout << "The linear problem has no solution." << endl;
		}

		cout << endl;
		
		/*
			Minimization problem
		*/
		objectiveFunction <<	3,
								4;

    // constraints = [A B], where A x >= B
		constraints <<		
      2,	1,	8,
			1,	2,	13,
		  1,	5,	16;
		
		solver2 = new SimplexSolver(SIMPLEX_MINIMIZE, objectiveFunction, constraints);

		if (solver2->hasSolution()) {
			cout << "The minimum is: " << solver2->getOptimum() << endl;
			cout << "The solution is: " << solver2->getSolution().transpose() << endl;
		} else {
			cout << "The linear problem has no solution." << endl;
		}
    VectorXd x = solver2->getSolution();
    cout<<(constraints.block(0,0,constraints.rows(),x.rows()) * x).transpose() <<
      " >?= "<<constraints.col(constraints.cols()-1).transpose()<<endl;
	}
	catch (exception *ex) {
    std::cout << ex->what() << std::endl;
	}
	
	delete solver1;
	delete solver2;
	
	return 0;
}

