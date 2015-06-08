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

// ---------------------------------------------------------
//
//  peng2004_surface_intergal_gradient_mex.cpp
//
// PENG2004_SURFACE_INTEGRAL_GRADIENT_MEX
// [int,grad_int] = peng2004_surface_integral_gradient_mex(X,V0,F0)
//
// Given a mesh (V0,F0), calculates the integral and gradient of the
// integral that defines the p-distance (p=3) from a set a points X
// to (V0,F0).
//
//  To compile:
//  mex peng2004_surface_integral_gradient_mex.cpp -lstdc++ -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include /Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/lib/libigl.a -I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
//
//  or rather:
//
//      mex LDFLAGS="\$LDFLAGS -framework Accelerate"  ...
//        peng2004_surface_integral_gradient_mex.cpp ...
//        -I/usr/local/igl/libigl/include ...
//        -I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported ...
//        /usr/local/igl/libigl/lib/libigl.a  -o peng2004_surface_integral_gradient_mex.cpp

#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
//#include <igl/linprog.h>
//#include <igl/mosek_linprog.h>
#include <igl/read_triangle_mesh.h>
#include <igl/circulation.h>
#include <igl/centroid.h>
#include <igl/matlab_format.h>
#include <igl/edges.h>
#include <igl/unique_edge_map.h>
#include <igl/unique.h>
#include <igl/intersect.h>
#include <igl/list_to_matrix.h>
#include <igl/matrix_to_list.h>
//#include <igl/Viewer/viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>
#include <algorithm>
#include <math.h>
#include <stdio.h>

#include <igl/write_triangle_mesh.h>

// Mex
#include <igl/matlab/MexStream.h>

#define MATLAB_LINK

using namespace std;
using namespace Eigen;
using namespace igl;

extern void _main();

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

void peng2004_canonical_wedge_integral_gradient(MatrixXd u,MatrixXd v, MatrixXd w,double beta,MatrixXd* int_1,MatrixXd* grad_int_1){

    int N = ((*int_1).rows());
    VectorXd q(N);
    double cotangent = 1.0/tan(beta/2);
    double sin_beta = sin(beta);
    double cos_beta = cos(beta);
    // norms of every vertex being moved
    for (int k=0; k<N; k++){
        q(k) = sqrt(u(k,0)*u(k,0)+v(k,0)*v(k,0)+w(k,0)*w(k,0));
        // Calculate integrals
        (*int_1)(k,0) = (2/w(k,0))*(atan2(w(k,0),(q(k,0)-u(k,0))*cotangent-v(k,0)));
        // Calculate gradient of the integral
        (*grad_int_1)(k,0) = sin_beta/(q(k,0)*(q(k,0)-(u(k,0)*cos_beta+v(k,0)*sin_beta)));
        (*grad_int_1)(k,1) = 1.0/(q(k,0)*(q(k,0)-u(k,0))) - cos_beta/(q(k,0)*(q(k,0)-(u(k,0)*cos_beta+v(k,0)*sin_beta)));
        (*grad_int_1)(k,2) = -(*int_1)(k,0)/w(k,0) +((u(k,0)*u(k,0)+v(k,0)*v(k,0) - q(k,0)*u(k,0))*sin_beta-q(k,0)*v(k,0)*(1.0-cos_beta))/(q(k,0)*w(k,0)*(q(k,0)-u(k,0))*(q(k,0)-(u(k,0)*cos_beta+v(k,0)*sin_beta)));
    }

    return;
}

void peng2004_wedge_integral_gradient(MatrixXd X,Vector3d P,Vector3d v1,Vector3d v2,MatrixXd* int_1,MatrixXd* grad_int_1){

    // first rotation: to make v1(2) = 0;
    double theta_z = atan2(v1(2),v1(1));
    MatrixXd R_theta_z(3,3);
    R_theta_z << 1.0, 0.0, 0.0,
                 0.0, cos(-theta_z), -sin(-theta_z),
                 0.0, sin(-theta_z), cos(-theta_z);
    // translate and rotate moving points
    X.col(0) = X.col(0).array()-P(0);
    X.col(1) = X.col(1).array()-P(1);
    X.col(2) = X.col(2).array()-P(2);
    X = (R_theta_z*X.transpose()).transpose();
    // rotate the wedge vectors
    v1 = R_theta_z*v1;
    v2 = R_theta_z*v2;
    // second rotation to make v1(1)=0
    double theta_y = atan2(v1(1),v1(0));
    MatrixXd R_theta_y(3,3);
    R_theta_y << cos(-theta_y), -sin(-theta_y), 0.0,
                 sin(-theta_y), cos(-theta_y), 0.0,
                 0.0, 0.0, 1.0;
    X = (R_theta_y*X.transpose()).transpose();
    // rotate the wedge vectors
    v1 = R_theta_y*v1;
    v2 = R_theta_y*v2;
    //  third rotation to make v2(2) = 0
    double theta_3 = atan2(v2(2),v2(1));
    MatrixXd R_theta_3(3,3);
    R_theta_3 << 1.0, 0.0, 0.0,
                 0.0, cos(-theta_3), -sin(-theta_3),
                 0.0, sin(-theta_3), cos(-theta_3);
    X = (R_theta_3*X.transpose()).transpose();
    v1 = R_theta_3*v1;
    v2 = R_theta_3*v2;
    // define beta angle - to belong to [0,2*pi]
    double beta = atan2(v2(1),v2(0));
    peng2004_canonical_wedge_integral_gradient(X.col(0),X.col(1),X.col(2),beta,int_1,grad_int_1);
    // back to original coordinate system
    MatrixXd R_theta_3_inv(3,3);
    R_theta_3_inv << 1.0, 0.0, 0.0,
                 0.0, cos(theta_3), -sin(theta_3),
                 0.0, sin(theta_3), cos(theta_3);
    *grad_int_1 = (R_theta_3_inv*(*grad_int_1).transpose()).transpose();
    MatrixXd R_theta_y_inv(3,3);
    R_theta_y_inv << cos(theta_y), -sin(theta_y), 0.0,
                 sin(theta_y), cos(theta_y), 0.0,
                 0.0, 0.0, 1.0;
    *grad_int_1 = (R_theta_y_inv*(*grad_int_1).transpose()).transpose();
    MatrixXd R_theta_z_inv(3,3);
    R_theta_z_inv << 1.0, 0.0, 0.0,
                 0.0, cos(theta_z), -sin(theta_z),
                 0.0, sin(theta_z), cos(theta_z);
    *grad_int_1 = (R_theta_z_inv*(*grad_int_1).transpose()).transpose();

    return;

}

void spatial_hash_faces(){


    return;

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

    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
                          "peng2004_surface_intergal_gradient_mex requires three input arguments.");
    }
    else if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
                          "peng2004_surface_intergal_gradient_mex generates two output arguments.");
    }

    const mxArray* X_mx = prhs[0];
    const mxArray* V0_mx = prhs[1];
    const mxArray* F0_mx  = prhs[2];

    vector<vector<double> > X_ = readMatrix<double>(X_mx);
    vector<vector<double> > V0_ = readMatrix<double>(V0_mx);
    vector<vector<double> > F0_  = readMatrix<double>(F0_mx);

    MatrixXd X,V0;
    MatrixXi F0;

    // Convert input to Eigen format
    if(V0_.size() > 0)
    {
        if(!list_to_matrix(V0_,V0))
        {
            return;
        }
        if(!list_to_matrix(X_,X))
        {
            return;
        }
        polygon_mesh_to_triangle_mesh(F0_,F0);
    }

    for (int i=0; i<F0.rows(); i++){
        for (int j=0; j<F0.cols(); j++){
            F0(i,j) = F0(i,j)-1;
        }
    }

    // This should be replaced by a spatial hash of face indices:
    // For each vertex in X, we would find to which cell it belongs,
    // collect all the faces in this cell and integrate over this cell

    // For each vertex compute a min_d value and gather faces that
    // are closer than (2*min_d+0.1)
    double min_d[X.rows()];
    vector<vector<int> > faces;
    faces.resize(X.rows());
    for (int i=0; i<X.rows();i++){
        min_d[i] = 1.0e+6;
        for (int j=0; j<F0.rows(); j++){
            // distance to the first vertex is the reference
            Vector3d v1(V0(F0(j,0),0),V0(F0(j,0),1),V0(F0(j,0),2));
            min_d[i] = fmin((v1-X.row(i).transpose()).norm(),min_d[i]);
        }
        for (int j=0; j<F0.rows(); j++){
            Vector3d v1(V0(F0(j,0),0),V0(F0(j,0),1),V0(F0(j,0),2));
            if ((v1-X.row(i).transpose()).norm()<2*min_d[i]+0.05){
                faces[i].push_back(j);
            }
        }

    }

//    double min_d = 1.0e+6;
//    for (int j=0; j<F0.rows(); j++){
//        Vector3d v1(V0(F0(j,0),0),V0(F0(j,0),1),V0(F0(j,0),2));
//        min_d = fmin((v1-X.row(0).transpose()).norm(),min_d);
//    }

    MatrixXd int_all = MatrixXd::Zero(X.rows(),1);
    MatrixXd grad_int_all = MatrixXd::Zero(X.rows(),3);

    // loop over points in X, loop over faces in its neighborhood
    for (int i=0; i<X.rows();i++){

        for (unsigned k=0; k<faces[i].size(); k++){

            int j = faces[i][k];
//            mexPrintf("i=%d, k=%d, j=%d \n", i,k,j);

            // Triangle vertices
            Vector3d v1(V0(F0(j,0),0),V0(F0(j,0),1),V0(F0(j,0),2));
            Vector3d v2(V0(F0(j,1),0),V0(F0(j,1),1),V0(F0(j,1),2));
            Vector3d v3(V0(F0(j,2),0),V0(F0(j,2),1),V0(F0(j,2),2));

            // First wedge
            Vector3d P_1 = v1;
            Vector3d e1_1 = v2-v1;
            Vector3d e2_1 = v3-v1;

            // Second wedge
            Vector3d P_2 = v3;
            Vector3d e1_2 = v2-v3;
            Vector3d e2_2 = v3-v1;

            // Third wedge
            Vector3d P_3 = v2;
            Vector3d e1_3 = v2-v3;
            Vector3d e2_3 = v2-v1;

            // Output for fisrt wedge
            MatrixXd int_1 = MatrixXd::Zero(1,1);
            MatrixXd grad_int_1(1,3);
            peng2004_wedge_integral_gradient(X.row(i),P_1,e1_1,e2_1,&int_1,&grad_int_1);

            // Output for second wedge
            MatrixXd int_2 = MatrixXd::Zero(1,1);
            MatrixXd grad_int_2(1,3);
            peng2004_wedge_integral_gradient(X.row(i),P_2,e1_2,e2_2,&int_2,&grad_int_2);

            // Output for third wedge
            MatrixXd int_3 = MatrixXd::Zero(1,1);
            MatrixXd grad_int_3(1,3);
            peng2004_wedge_integral_gradient(X.row(i),P_3,e1_3,e2_3,&int_3,&grad_int_3);

            int_all.row(i) = int_all.row(i) + (int_1-int_2+int_3);
            grad_int_all.row(i) = grad_int_all.row(i) + (grad_int_1-grad_int_2+grad_int_3);


        }

    }

//    MatrixXd int_all = MatrixXd::Zero(X.rows(),1);
//    MatrixXd grad_int_all = MatrixXd::Zero(X.rows(),3);


//    // comment this to have a code to loop over vertices
//    for (int j=0; j<F0.rows(); j++){

//        // Triangle vertices
//        Vector3d v1(V0(F0(j,0),0),V0(F0(j,0),1),V0(F0(j,0),2));
//        Vector3d v2(V0(F0(j,1),0),V0(F0(j,1),1),V0(F0(j,1),2));
//        Vector3d v3(V0(F0(j,2),0),V0(F0(j,2),1),V0(F0(j,2),2));

//        // add here a parameter max_d such that the calculation is done
//        // only if the vertex is closer than max_d. This indeed speeds it up
////        if ((v1-X.row(0).transpose()).norm()<5*min_d){

//            // First wedge
//            Vector3d P_1 = v1;
//            Vector3d e1_1 = v2-v1;
//            Vector3d e2_1 = v3-v1;

//            // Second wedge
//            Vector3d P_2 = v3;
//            Vector3d e1_2 = v2-v3;
//            Vector3d e2_2 = v3-v1;

//            // Third wedge
//            Vector3d P_3 = v2;
//            Vector3d e1_3 = v2-v3;
//            Vector3d e2_3 = v2-v1;

//            // Output for fisrt wedge
//            MatrixXd int_1 = MatrixXd::Zero(X.rows(),1);
//            MatrixXd grad_int_1(X.rows(),3);
//            peng2004_wedge_integral_gradient(X,P_1,e1_1,e2_1,&int_1,&grad_int_1);

//            // Output for second wedge
//            MatrixXd int_2 = MatrixXd::Zero(X.rows(),1);
//            MatrixXd grad_int_2(X.rows(),3);
//            peng2004_wedge_integral_gradient(X,P_2,e1_2,e2_2,&int_2,&grad_int_2);

//            // Output for third wedge
//            MatrixXd int_3 = MatrixXd::Zero(X.rows(),1);
//            MatrixXd grad_int_3(X.rows(),3);
//            peng2004_wedge_integral_gradient(X,P_3,e1_3,e2_3,&int_3,&grad_int_3);

//            int_all = int_all + (int_1-int_2+int_3);
//            grad_int_all = grad_int_all + (grad_int_1-grad_int_2+grad_int_3);

////        }

//    }

    vector<vector<double> > int_all_;
    vector<vector<double> > grad_int_all_;
    matrix_to_list(int_all,int_all_);
    matrix_to_list(grad_int_all,grad_int_all_);

    // Output to Matlab
    plhs[0] = writeMatrix(int_all_);
    plhs[1] = writeMatrix(grad_int_all_);

    // Obtaining different values for different runs (check it)

    // delete stuff
    X_.clear();
    V0_.clear();
    F0_.clear();

    // Restore the std stream buffer Important!
    std::cout.rdbuf(outbuf);
    return;

    }
