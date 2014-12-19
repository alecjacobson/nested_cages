
#include <iostream>
#include <math.h>
#include "mex.h"
#include <vector>
#include <fstream>

#include <VelocityFilter.h>

#define MATLAB_LINK

using namespace std;
using namespace Eigen;

extern void _main();

void parse_rhs(
  const int nrhs,
  const mxArray *prhs[],
  Eigen::MatrixXd & V)
{
  using namespace std;
  // set number of mesh vertices
  const int n = mxGetM(prhs[0]);
  // set vertex position pointers
  double * Vp = mxGetPr(prhs[0]);
  const int dim = mxGetN(prhs[0]);
  // resize output to transpose
  V.resize(n,dim);
  copy(Vp,Vp+n*dim,&V.data()[0]);
}

void parse_rhs(
  const int nrhs,
  const mxArray *prhs[],
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F)
{
  using namespace std;
  if(nrhs < 2)
  {
    mexErrMsgTxt("nrhs < 2");
  }

  parse_rhs(nrhs,prhs,V);

  const int dim = V.cols();
  if(dim != 3 && dim != 2)
  {
    mexErrMsgTxt("Mesh vertex list must be #V by 2 or 3 list of vertex positions");
  }
  if(dim != (int)mxGetN(prhs[1]))
  {
   mexErrMsgTxt("Mesh facet size must equal dimension");
  }

  // set number of faces
  const int m = mxGetM(prhs[1]);
  // set face index list pointer
  double * Fp = mxGetPr(prhs[1]);

  if((int)mxGetN(prhs[1]) != dim)
  {
    mexErrMsgTxt("Origin list dimension must match vertex list dimension");
  }

  // resize output to transpose
  F.resize(m,dim);
  // Q: Is this doing a cast?
  // A: Yes.
  copy(Fp,Fp+m*dim,F.data());
  // http://stackoverflow.com/a/4461466/148668
  transform(F.data(),F.data()+m*dim,F.data(),
    bind2nd(std::plus<double>(),-1.0));
}

void parse_rhs(
  const int nrhs,
  const mxArray *prhs[],
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & U,
  Eigen::MatrixXi & G)
{
  parse_rhs(nrhs,prhs,V,F);
  parse_rhs(nrhs-2,prhs+2,U,G);
}

void parse_rhs(
  const int nrhs,
  const mxArray *prhs[],
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F)
{
  parse_rhs(nrhs,prhs,P);
  parse_rhs(nrhs-1,prhs+1,V,F);
}

void parse_rhs(
  const int nrhs,
  const mxArray *prhs[],
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
        int *extra_number)
{
  parse_rhs(nrhs,prhs,P);
  parse_rhs(nrhs-1,prhs+1,V,F);
  *extra_number = *mxGetPr(prhs[3]);
}


void mexFunction(
		 int          nlhs,
		 mxArray      *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{

  /* Check for proper number of arguments */

  if (nrhs != 4) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "velocity_filter requires four input arguments.");
  } 
  else if (nlhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "velocity_filter generates two one argument.");
  }

  MatrixXd V0, V1;
  MatrixXi F0;
  int numinfinite;
  parse_rhs(nrhs,prhs,V0,V1,F0,&numinfinite);
  Matrix3Xi F0_t(3,F0.rows());
  for (int k=0; k<F0.rows();k++){
      F0_t(0,k) = F0(k,0);
      F0_t(1,k) = F0(k,1);
      F0_t(2,k) = F0(k,2);
  }

  VectorXd qstart(3*V0.rows());
  VectorXd qend(3*V0.rows());
  for (int k=0; k<V0.rows();k++){
      qstart(3*k) = V0(k,0);
      qstart(3*k+1) = V0(k,1);
      qstart(3*k+2) = V0(k,2);
  }
  for (int k=0; k<V0.rows();k++){
      qend(3*k) = V1(k,0);
      qend(3*k+1) = V1(k,1);
      qend(3*k+2) = V1(k,2);
  }
  VectorXd invmasses(3*V0.rows());
  for (int k=0; k<3*numinfinite; k++){
      invmasses(k) = 0.0;
  }
  for (int k=3*numinfinite; k<3*V0.rows(); k++){
      invmasses(k) = 1.0;
  }

  int sim_status = VelocityFilter::velocityFilter(qstart, qend, F0_t, invmasses, 1e-4, 1e-6);
  mexPrintf("simulation_status = %d \n", sim_status);

  for (int k=0; k<V0.rows(); k++){
      V1(k,0) = qend(3*k);
      V1(k,1) = qend(3*k+1);
      V1(k,2) = qend(3*k+2);
  }

  plhs[0] = mxCreateDoubleMatrix(V1.rows(),V1.cols(), mxREAL);
  double * Vp = mxGetPr(plhs[0]);
  copy(&V1.data()[0],&V1.data()[0]+V1.size(),Vp);

  return;
}
