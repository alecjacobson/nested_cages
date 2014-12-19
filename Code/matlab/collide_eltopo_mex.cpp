// ---------------------------------------------------------
//
//  collide_eltopo_mex.cpp
//  Leonardo Sacht, 2014.
//
// COLLIDE_ELTOPO_MEX
// collide_eltopo_mex(V0,F0,V1,N)
//
// Takes 2 meshes (V0,F0) and (V1,F0) (with the same connectivity)
// and a number N that defines the first vertices that will have
// infinite mass.
//
//  To compile:
//  mex collide_eltopo_mex.cpp -I../common -I../eltopo3d -I../talpa -I../talpa/drivers -I../common/tunicate -llapack -lblas -lstdc++ /Users/Leo/PHD_Work/Cage_Generation_2013/code/starlab-mcfskel/core/external/cholmod-4.0.0/lib/osx64/libcholmod.a /Users/Leo/PHD_Work/Cage_Generation_2013/code/starlab-mcfskel/core/external/cholmod-4.0.0/lib/osx64/libamd.a /Users/Leo/PHD_Work/Cage_Generation_2013/code/starlab-mcfskel/core/external/cholmod-4.0.0/lib/osx64/libcolamd.a libeltopo_release.a -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include/igl /opt/local/lib/gcc47/libgomp.a /Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/lib/libigl.a /opt/local/lib/libSuiteSparse.dylib ../talpa/obj/bfstream.o ../talpa/obj/iomesh.o
//
//  or rather:
//
//      mex LDFLAGS="\$LDFLAGS -framework Accelerate"  ...
//        collide_eltopo_mex.cpp -I../eltopo/common -I../eltopo/eltopo3d ...
//        -I../eltopo/talpa -I../eltopo/talpa/drivers ...
//        -I../eltopo/common/tunicate -L../eltopo/eltopo3d/ ...
//        -leltopo_release -o collide_eltopo_mex

// std
#include <cstdio>
#include <fenv.h>
#include <fstream>
#include <vector>
#include <queue>
#include <cfloat>

// common
#include <array2.h>
#include <ccd_wrapper.h>
#include <collisionqueries.h>
#include <expansion.h>
#include <marching_tiles_hires.h>
#include <util.h>
#include <vec.h>
#include <wallclocktime.h>

// el topo
#include <collisionpipeline.h>
#include <eltopo.h>
#include <iomesh.h>
#include <meshrenderer.h>
#include <runstats.h>
#include <surftrack.h>
#include <trianglequality.h>

// talpa
#include <framestepper.h>
#include <meancurvature.h>
#include <meshdriver.h>
#include <scriptinit.h>
#include <simulation.h>

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include <iostream>
#include <math.h>
#include "mex.h"
#include <vector>
#include <fstream>

#define MATLAB_LINK

using namespace std;

ofstream* os;
vector<vector<double> > V_eltopo;

extern void _main();

int vertices_fine;
double eps_prox;
double tol_dt;

template<class T>
vector<vector<T> > readMatrix(const mxArray* mat)
{
    T* ptr = mxGetPr(mat);

    int m = mxGetM(mat);
    int n = mxGetN(mat);

    vector<vector<T> >	V;
    V.resize(m);
        for(int j=0;j<n;j++)
            for(int i=0;i<m;++i)
                V[i].push_back(*(ptr++));

    return V;

}

template<class T>
mxArray* writeMatrix(vector<vector<T> > mat)
{
    mxArray* ret = 0;

    if (mat.size() == 0)
    {
        ret =	mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    }
    else
    {
        int M = mat.size();
        int N = mat[0].size();

        ret = mxCreateNumericMatrix(M, N, mxDOUBLE_CLASS, mxREAL);
    double* pointer = mxGetPr(ret);

    /* Copy data into the mxArray */
    int count = 0;
    for ( int j = 0; j < N; ++j )
    {
        for (int i = 0; i < M; ++i)
        {
            pointer[count++] = mat[i][j];
        }
    }
  }

  return ret;
}

// ---------------------------------------------------------
// Global interface declarations
// ---------------------------------------------------------

extern "C" {
    void exactinit();    // in predicates.c
}

void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{

  V_eltopo.clear();

//  cout << "Entered mexFunction" << endl;

  /* Check for proper number of arguments */

  if (nrhs != 6) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
            "collide_eltopo_mex six input arguments.");
  }
  else if (nlhs != 2) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "collide_eltopo_mex generates two output arguments.");
  }

  const mxArray* V0_mx = prhs[0];
  const mxArray* F_mx = prhs[1];
  const mxArray* V1_mx  = prhs[2];
  // get number of vertices of the fine mesh
  vertices_fine = *mxGetPr(prhs[3]);
  eps_prox = *mxGetPr(prhs[4]);
  tol_dt = *mxGetPr(prhs[5]);

  os = new ofstream("/Users/Leo/log.txt");

  vector<vector<double> > V0_ = readMatrix<double>(V0_mx);
  vector<vector<double> > V1_ = readMatrix<double>(V1_mx);
  vector<vector<double> > F  = readMatrix<double>(F_mx );

  // Mesh at time 0
  // get number of vertices
  int num_vertices = mxGetM(V0_mx);
  // get vertex positions
  double V0[3*num_vertices];
  for (int i = 0; i<num_vertices; i++){
      V0[3*i] = V0_[i][0];
      V0[3*i+1] = V0_[i][1];
      V0[3*i+2] = V0_[i][2];
  }
  // get number of triangles
  int num_triangles = mxGetM(F_mx);
  // get face indices
  int F0[3*num_triangles];
  for (int i = 0; i<num_triangles; i++){
      F0[3*i] = F[i][0]-1;
      F0[3*i+1] = F[i][1]-1;
      F0[3*i+2] = F[i][2]-1;
  }
  // set vertex masses (=infty for the first vertices and =eps for the last ones)
  double masses[num_vertices];
  for (int i=0; i<vertices_fine; i++){
      masses[i] = std::numeric_limits<double>::infinity();
  }
  for (int i=vertices_fine; i<num_vertices; i++){
      masses[i] = 1.0;
  }
  // encapsulate all data into an ElTopoMesh
  ElTopoMesh eltopo_time0;
  eltopo_time0.num_vertices = num_vertices;
  eltopo_time0.vertex_locations = V0;
  eltopo_time0.num_triangles = num_triangles;
  eltopo_time0.triangles = F0;
  eltopo_time0.vertex_masses = masses;

  // Mesh at time 1
  double V1[3*num_vertices];
  // get vertex positions
  for (int i = 0; i<num_vertices; i++){
      V1[3*i] = V1_[i][0];
      V1[3*i+1] = V1_[i][1];
      V1[3*i+2] = V1_[i][2];
  }

  // Set general parameters
  ElTopoGeneralOptions sim_general_options;
  // do not print stuff to the console
  sim_general_options.m_verbose = 0;
  // do avoid self-intersections
  sim_general_options.m_collision_safety = 1;
  sim_general_options.m_proximity_epsilon = eps_prox;

  // Set Simulation parameters
  ElTopoIntegrationOptions sim_integration_options;
  sim_integration_options.m_friction_coefficient = 0.0;
  sim_integration_options.m_dt = 1.0;

  // run simulation (cut time step if it does not satisfy constraints)
  double* V_final;
  double out_dt = 0.0;
  double rest_dt = 1.0;
  int attempts = 0;
  while (rest_dt > tol_dt){
      el_topo_integrate(&eltopo_time0, V1, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
      eltopo_time0.vertex_locations = V_final;
      mexPrintf("out_dt = %.4f\n", out_dt);
      rest_dt = (1-out_dt)*rest_dt;
      mexPrintf("rest_dt = %.4f\n", rest_dt);
      attempts = attempts+1;

      if (out_dt<tol_dt){
        // return something as output
        plhs[0] = writeMatrix(V0_);
        // it didn't find a time step, back to Matlab prompt
        plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
        double *out_ptr;
        out_ptr = mxGetPr(plhs[1]);
        *out_ptr= rest_dt;
        return;
      }
  }

  // get vertex positions
  for (int i = 0; i<num_vertices; i++){
      vector<double> temp(3);
      temp[0] = V_final[3*i];
      temp[1] = V_final[3*i+1];
      temp[2] = V_final[3*i+2];

      V_eltopo.push_back(temp);
  }

  os->close();

  plhs[0] = writeMatrix(V_eltopo);
  plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
  double *out_ptr;
  out_ptr = mxGetPr(plhs[1]);
  *out_ptr= rest_dt;

  V_eltopo.clear();

  return;
}

