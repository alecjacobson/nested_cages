#include <iostream>
#include <math.h>
#include "mex.h"
#include <vector>
#include <fstream>

#include <VelocityFilter.h>
#include <Distance.h>
#include <igl/matlab/MexStream.h>

#define MATLAB_LINK

using namespace std;
using namespace Eigen;

ofstream* os;

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
  igl::MexStream mout;
  std::streambuf *outbuf = cout.rdbuf(&mout);

  /* Check for proper number of arguments */

  if (nrhs != 7) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
            "testNewSequence requires seven input arguments.");
  }
  else if (nlhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "testNewSequence generates one output argument.");
  }

  os = new ofstream("./testNewSequence_mex-log.txt");

  MatrixXd V_coarse,Pall;
  MatrixXi F_coarse,F_fine;
  parse_rhs(nrhs,prhs,V_coarse,F_coarse,Pall,F_fine);

  int numfineverts = *mxGetPr(prhs[4]);
  double outerRadius = *mxGetPr(prhs[5]);
  double innerRadius = *mxGetPr(prhs[6]);

  // Adapting Matlab's coarse mesh to Etienne's code
  VectorXd qcoarse(3*V_coarse.rows());
  for (int k=0; k<V_coarse.rows();k++){
      qcoarse(3*k) = V_coarse(k,0);
      qcoarse(3*k+1) = V_coarse(k,1);
      qcoarse(3*k+2) = V_coarse(k,2);
  }
  int numcoarseverts = qcoarse.size()/3;
  Matrix3Xi fcoarse(3,F_coarse.rows());
  for (int k=0; k<F_coarse.rows();k++){
      fcoarse(0,k) = F_coarse(k,0);
      fcoarse(1,k) = F_coarse(k,1);
      fcoarse(2,k) = F_coarse(k,2);
  }
  mexPrintf("Loaded coarse mesh with %d vertices and %d faces \n", numcoarseverts, fcoarse.cols());

  // Adapting Matlab's (initial) fine mesh to Etienne's code
  VectorXd qfine(3*numfineverts);
  int numfinemeshes = Pall.rows()/numfineverts;
  mexPrintf("Number of fine meshes is %d \n", numfinemeshes);
  for (int k=0; k<numfineverts;k++){
      qfine(3*k) = Pall((numfinemeshes-1)*numfineverts+k,0);
      qfine(3*k+1) = Pall((numfinemeshes-1)*numfineverts+k,1);
      qfine(3*k+2) = Pall((numfinemeshes-1)*numfineverts+k,2);
  }
  Matrix3Xi ffine(3,F_fine.rows());
  for (int k=0; k<F_fine.rows();k++){
      ffine(0,k) = F_fine(k,0);
      ffine(1,k) = F_fine(k,1);
      ffine(2,k) = F_fine(k,2);
  }
  mexPrintf("Loaded fine mesh with %d vertices and %d faces \n", qfine.size()/3, ffine.cols());

  // Now copy and paste from Etienne's code
  VectorXd q(qcoarse.size() + qfine.size());
  Matrix3Xi f;
  f.resize(3, fcoarse.cols() + ffine.cols());
  for(int i=0; i<(int)qcoarse.size(); i++)
      q[i] = qcoarse[i];
  for(int i=0; i<qfine.size(); i++)
      q[qcoarse.size()+i] = qfine[i];
  for(int i=0; i<(int)fcoarse.cols(); i++)
      f.col(i) = fcoarse.col(i);
  for(int i=0; i<(int)ffine.cols(); i++)
      for(int j=0; j<3; j++)
          f.coeffRef(j, fcoarse.cols()+i) = ffine.coeff(j,i) + numcoarseverts;

  VectorXd invmasses(qcoarse.size() + qfine.size());
  for(int i=0; i<(int)qcoarse.size(); i++)
      invmasses[i] = 1.0;
  for(int i=0; i<(int)qfine.size(); i++)
      invmasses[qcoarse.size()+i] = 0;

  set<int> fixedVerts;
  for(int i=0; i<(int)qfine.size()/3; i++)
      fixedVerts.insert(qcoarse.size()/3 + i);

  double distance = Distance::meshSelfDistance(q, f, fixedVerts);
  while(distance < outerRadius)
  {
      mexPrintf("Distance is %.16f. Inflating ... \n", distance);
      VectorXd qnew = q;
      VelocityFilter::velocityFilter(q, qnew, f, invmasses, 2.0*distance, 0.5*distance);
      q = qnew;
      distance = Distance::meshSelfDistance(q, f, fixedVerts);
  }
  for(int i=0; i<(int)qcoarse.size(); i++)
      qcoarse[i] = q[i];

  int curmesh = 0;
  while(curmesh<numfinemeshes-1)
  {

      mexPrintf("\n \n \n Now processing Mesh %d -> Mesh %d. Total is %d \n", curmesh, curmesh+1, numfinemeshes-1);

      VectorXd qfine1(3*numfineverts);
      for (int k=0; k<numfineverts;k++){
          qfine1(3*k) = Pall((numfinemeshes-1-curmesh)*numfineverts+k,0);
          qfine1(3*k+1) = Pall((numfinemeshes-1-curmesh)*numfineverts+k,1);
          qfine1(3*k+2) = Pall((numfinemeshes-1-curmesh)*numfineverts+k,2);
      }
      Matrix3Xi ffine1 = ffine;
      VectorXd qfine2(3*numfineverts);
      for (int k=0; k<numfineverts;k++){
          qfine2(3*k) = Pall((numfinemeshes-1-(curmesh+1))*numfineverts+k,0);
          qfine2(3*k+1) = Pall((numfinemeshes-1-(curmesh+1))*numfineverts+k,1);
          qfine2(3*k+2) = Pall((numfinemeshes-1-(curmesh+1))*numfineverts+k,2);
      }
      Matrix3Xi ffine2 = ffine;

      mexPrintf("Loaded fine mesh with %d vertices and %d faces \n", qfine1.size()/3, ffine1.cols());

      VectorXd q1(qcoarse.size() + qfine1.size());
      VectorXd q2(qcoarse.size() + qfine2.size());
      Matrix3Xi f1;
      f1.resize(3, fcoarse.cols() + ffine1.cols());
      Matrix3Xi f2;
      f2.resize(3, fcoarse.cols() + ffine2.cols());
      for(int i=0; i<(int)qcoarse.size(); i++)
      {
          q1[i] = qcoarse[i];
          q2[i] = qcoarse[i];
      }
      for(int i=0; i<(int)fcoarse.cols(); i++)
      {
          f1.col(i) = fcoarse.col(i);
          f2.col(i) = fcoarse.col(i);
      }
      for(int i=0; i<qfine1.size(); i++)
          q1[qcoarse.size() + i] = qfine1[i];
      for(int i=0; i<qfine2.size(); i++)
          q2[qcoarse.size() + i] = qfine2[i];
      for(int i=0; i<(int)ffine1.cols(); i++)
      {
          for(int j=0; j<3; j++)
              f1.coeffRef(j, fcoarse.cols() + i) = ffine1.coeff(j, i) + numcoarseverts;
      }
      for(int i=0; i<(int)ffine2.cols(); i++)
      {
          for(int j=0; j<3; j++)
              f2.coeffRef(j, fcoarse.cols() + i) = ffine2.coeff(j, i) + numcoarseverts;
      }
      VectorXd invmasses(qcoarse.size() + qfine1.size());
      for(int i=0; i<(int)qcoarse.size(); i++)
          invmasses[i] = 1.0;
      for(int i=0; i<(int)qfine1.size(); i++)
          invmasses[qcoarse.size()+i] = 0;

      int ret = VelocityFilter::velocityFilter(q1, q2, f1, invmasses, outerRadius, innerRadius);
      if(ret < 0)
      {
          mexPrintf("Velocity filter failed! Error code: %d", ret);
          return;
      }
      for(int i=0; i<(int)qcoarse.size(); i++)
          qcoarse[i] = q2[i];

//      stringstream ss;
//      ss << "out-iter-" << curmesh << ".obj";

//      writeMesh(ss.str().c_str(), qcoarse, fcoarse);
      curmesh++;
  }
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);

  // Output final coarse mesh
  MatrixXd V1(V_coarse.rows(),3);
  for (int k=0; k<V_coarse.rows(); k++){
      V1(k,0) = qcoarse(3*k);
      V1(k,1) = qcoarse(3*k+1);
      V1(k,2) = qcoarse(3*k+2);
  }

  plhs[0] = mxCreateDoubleMatrix(V1.rows(),V1.cols(), mxREAL);
  double * Vp = mxGetPr(plhs[0]);
  copy(&V1.data()[0],&V1.data()[0]+V1.size(),Vp);

  return;

}
