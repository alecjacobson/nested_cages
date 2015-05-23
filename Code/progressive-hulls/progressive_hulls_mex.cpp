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
//  progressive_hulls_mex.cpp
//
// PROGRESSIVE_HULLS_MEX
// [hulls_V,hulls_F] = progressive_hulls_mex(V0,F0,levels,w)
//
// Given an input mesh (V0,F0), computes the Progressive Hulls
// with resolution specified by the vector levels. The Progressive
// Hulls are descibed in "Silhouette Clipping", Sander et al. [2000].
//
//  To compile:
//  mex progressive_hulls_mex.cpp -lstdc++ -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include /Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/lib/libigl.a -I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
//
//  or rather:
//
//      mex LDFLAGS="\$LDFLAGS -framework Accelerate"  ...
//        collide_eltopo_mex.cpp -I../eltopo/common -I../eltopo/eltopo3d ...
//        -I/usr/local/igl/libigl/include ...
//        -I../eltopo/talpa -I../eltopo/talpa/drivers ...
//        -I../eltopo/common/tunicate -L../eltopo/eltopo3d/ ...
//        -leltopo_release -o collide_eltopo_mex

#include "eiquadprog.hpp"
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

#include <igl/write_triangle_mesh.h>

// Mex
#include <igl/matlab/MexStream.h>

#define MATLAB_LINK

using namespace std;

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

// For debuggin'
int at(
  Eigen::MatrixXi & M,
  const int i,
  const int j)
{
  return M(i,j);
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

    using namespace Eigen;
    using namespace igl;

    /* Check for proper number of arguments */

    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
                          "progressive_hulls_mex requires four input arguments.");
    }
    else if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
                          "progressive_hulls_mex generates two output arguments.");
    }

    const mxArray* V0_mx = prhs[0];
    const mxArray* F0_mx = prhs[1];
    const mxArray* levels_mx  = prhs[2];

    double w = mxGetScalar(prhs[3]);;

    vector<vector<double> > V0_ = readMatrix<double>(V0_mx);
    vector<vector<double> > F0_  = readMatrix<double>(F0_mx);
    vector<vector<double> > levels_ = readMatrix<double>(levels_mx);

    MatrixXd V,OV;
    MatrixXi F,OF;

    // Convert input to Eigen format
    if(V0_.size() > 0)
    {
        if(!list_to_matrix(V0_,OV))
        {
            return;
        }
        polygon_mesh_to_triangle_mesh(F0_,OF);
    }

    for (int i=0; i<OF.rows(); i++){
        for (int j=0; j<OF.cols(); j++){
            OF(i,j) = OF(i,j)-1;
        }
    }


    // Here comes Alec's code
    // first produce individual prog hulls, then output all of them

    VectorXi EMAP;
    MatrixXi E,EF,EI;
    typedef std::set<std::pair<double,int> > PriorityQueue;
    PriorityQueue Q;
    std::vector<PriorityQueue::iterator > Qit;
    // If an edge were collapsed, we'd collapse it to these points:
    MatrixXd C;
    int num_collapsed;

    std::function<void(
                const int,
                const Eigen::MatrixXd &,
                const Eigen::MatrixXi &,
                const Eigen::MatrixXi &,
                const Eigen::VectorXi &,
                const Eigen::MatrixXi &,
                const Eigen::MatrixXi &,
                double &,
                RowVectorXd &)> cost_and_placement;
    const auto & shortest_edge_and_midpoint = [](
            const int e,
            const Eigen::MatrixXd & V,
            const Eigen::MatrixXi & /*F*/,
            const Eigen::MatrixXi & E,
            const Eigen::VectorXi & /*EMAP*/,
            const Eigen::MatrixXi & /*EF*/,
            const Eigen::MatrixXi & /*EI*/,
            double & cost,
            RowVectorXd & p)
    {
        cost = (V.row(E(e,0))-V.row(E(e,1))).norm();
        p = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
    };


    // Leo: This is the function to compute the cost of an edge collapse
    const auto & minimal_volume_and_outer_hull = [
            //&env
            &w
            ](
                const int e,
                const Eigen::MatrixXd & V,
                const Eigen::MatrixXi & F,
                const Eigen::MatrixXi & E,
                const Eigen::VectorXi & EMAP,
                const Eigen::MatrixXi & EF,
                const Eigen::MatrixXi & EI,
                double & cost,
                RowVectorXd & p)
    {
            assert(V.cols() == 3 && "V.cols() should be 3");
            // Gather list of unique face neighbors
            vector<int> Nall =  circulation(e, true,F,E,EMAP,EF,EI);
            vector<int> Nother= circulation(e,false,F,E,EMAP,EF,EI);
            Nall.insert(Nall.end(),Nother.begin(),Nother.end());
            vector<int> N;
            igl::unique(Nall,N);
            // Gather:
            //   A  #N by 3 normals scaled by area,
            //   D  #N determinants of matrix formed by points as columns
            //   B  #N point on plane dot normal
            MatrixXd A(N.size(),3);
            VectorXd D(N.size());
            VectorXd B(N.size());
            //cout<<"N=[";
            for(int i = 0;i<N.size();i++)
    {
        const int f = N[i];
        //cout<<(f+1)<<" ";
        const RowVector3d & v01 = V.row(F(f,1))-V.row(F(f,0));
        const RowVector3d & v20 = V.row(F(f,2))-V.row(F(f,0));
        A.row(i) = v01.cross(v20);
        B(i) = V.row(F(f,0)).dot(A.row(i));
        D(i) =
                (Matrix3d()<< V.row(F(f,0)), V.row(F(f,1)), V.row(F(f,2)))
                .finished().determinant();
    }
    //cout<<"];"<<endl;
    // linear objective
    Vector3d f = A.colwise().sum().transpose();
    VectorXd x;
    //bool success = linprog(f,-A,-B,MatrixXd(0,A.cols()),VectorXd(0,1),x);
    //VectorXd _;
    //bool success = mosek_linprog(f,A.sparseView(),B,_,_,_,env,x);
    //if(success)
    //{
    //  cost  = (1./6.)*(x.dot(f) - D.sum());
    //}
    bool success = false;
    {
        RowVectorXd mid = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
        MatrixXd G =  w*Matrix3d::Identity(3,3);
        VectorXd g0 = (1.-w)*f - w*mid.transpose();
        const int n = A.cols();
        success = solve_quadprog(
                    G,g0,
                    MatrixXd(n,0),VectorXd(0,1),
                    A.transpose(),-B,x);
        cost  = (1.-w)*(1./6.)*(x.dot(f) - D.sum()) +
                w*(x.transpose()-mid).squaredNorm() +
                w*(V.row(E(e,0))-V.row(E(e,1))).norm();
    }

    // A x >= B
    // A x - B >=0
    // This is annoyingly necessary. Seems the solver is letting some garbage
    // slip by.
    success = success && ((A*x-B).minCoeff()>-1e-10);
    if(success)
    {
        p = x.transpose();
        //assert(cost>=0 && "Cost should be positive");
    }else
    {
        cost = std::numeric_limits<double>::infinity();
        //VectorXi NM;
        //igl::list_to_matrix(N,NM);
        //cout<<matlab_format((NM.array()+1).eval(),"N")<<endl;
        //cout<<matlab_format(f,"f")<<endl;
        //cout<<matlab_format(A,"A")<<endl;
        //cout<<matlab_format(B,"B")<<endl;
        //exit(-1);
        p = RowVectorXd::Constant(1,3,std::nan("inf-cost"));
    }
};

  //cost_and_placement = shortest_edge_and_midpoint;
  cost_and_placement = minimal_volume_and_outer_hull;

  // Leo: Here it constructs the initial queue of edges
  const auto & reset = [&]()
  {
    F = OF;
    V = OV;
    edge_flaps(F,E,EMAP,EF,EI);
    Qit.resize(E.rows());

    C.resize(E.rows(),V.cols());
    VectorXd costs(E.rows());
    for(int e = 0;e<E.rows();e++)
    {
      double cost = e;
      RowVectorXd p(1,3);
      cost_and_placement(e,V,F,E,EMAP,EF,EI,cost,p);
      C.row(e) = p;
      Qit[e] = Q.insert(std::pair<double,int>(cost,e)).first;
    }
    num_collapsed = 0;
//    viewer.data.clear();
//    viewer.data.set_mesh(V,F);
//    viewer.data.set_face_based(true);
  };


  //  const auto &pre_draw = [&](igl::Viewer & viewer)->bool
    const auto &pre_draw = [&]()
    {
  //    if(viewer.core.is_animating && !Q.empty())
      if(!Q.empty())
      {
        bool something_collapsed = false;
        // collapse edge
  //      const int max_iter = std::ceil(0.01*Q.size());
        const int max_iter = 1; // here we collapse one edge at a time
        for(int j = 0;j<max_iter;j++)
        {
          if(Q.empty())
          {
            break;
          }
          std::pair<double,int> p = *(Q.begin());
          if(p.first == std::numeric_limits<double>::infinity())
          {
            break;
          }
          Q.erase(Q.begin());
          Qit[p.second] = Q.end();
          std::vector<int> N  = circulation(p.second, true,F,E,EMAP,EF,EI);
          std::vector<int> Nd = circulation(p.second,false,F,E,EMAP,EF,EI);
          N.insert(N.begin(),Nd.begin(),Nd.end());
          int e1,e2,f1,f2;
          const bool collapsed =
            collapse_edge(p.second,C.row(p.second),V,F,E,EMAP,EF,EI,e1,e2,f1,f2);
          if(collapsed)
          {
            something_collapsed = true;
            // Erase the two, other collapsed edges
            Q.erase(Qit[e1]);
            Qit[e1] = Q.end();
            Q.erase(Qit[e2]);
            Qit[e2] = Q.end();
            // update local neighbors
            // loop over original face neighbors
            for(auto n : N)
            {
              if(F(n,0) != IGL_COLLAPSE_EDGE_NULL ||
                 F(n,1) != IGL_COLLAPSE_EDGE_NULL ||
                 F(n,2) != IGL_COLLAPSE_EDGE_NULL)
              {
                for(int v = 0;v<3;v++)
                {
                  // get edge id
                  const int e = EMAP(v*F.rows()+n);
                  // erase old entry
                  Q.erase(Qit[e]);
                  // compute cost and potential placement
                  double cost;
                  RowVectorXd p;
                  cost_and_placement(e,V,F,E,EMAP,EF,EI,cost,p);
                  // Replace in queue
                  Qit[e] = Q.insert(std::pair<double,int>(cost,e)).first;
                  C.row(e) = p;
                }
              }
            }
            num_collapsed++;
          }else
          {
            // reinsert with infinite weight
            p.first = std::numeric_limits<double>::infinity();
            Qit[p.second] = Q.insert(p).first;
          }
        }
  //      if(something_collapsed)
  //      {
  //        viewer.data.clear();
  //        viewer.data.set_mesh(V,F);
  //        viewer.data.set_face_based(true);
  //      }
      }
  //    return false;
    };


      reset();
      int num_levels = mxGetN(levels_mx);
      mxArray *hulls_V_ptr, *hulls_F_ptr;
      hulls_V_ptr = mxCreateCellMatrix( num_levels+1, 1 );
      hulls_F_ptr = mxCreateCellMatrix( num_levels+1, 1 );
      // vector<vector>> format
      vector<vector<double> > Vf_;
      vector<vector<double> > Ff_;
      // output initial meshes as last level
      mxSetCell(hulls_V_ptr, num_levels, writeMatrix(V0_));
      mxSetCell(hulls_F_ptr, num_levels, writeMatrix(F0_));
      int cur_num_levels = 0;
      while (cur_num_levels<num_levels){

          int to_collapse = 0.5*(mxGetM(F0_mx)-levels_[0][num_levels-1-cur_num_levels]);
//          int to_collapse = 10;

          while (num_collapsed<to_collapse){
              pre_draw();
          }


        // And output to Matlab
        // (placeholder to test: output the input)
        MatrixXd Vf = V;
        MatrixXd Ff(F.rows(),F.cols());
        // Convert to double so Matlab can use it
        for (int i=0; i<F.rows(); i++){
            for (int j=0; j<F.cols(); j++){
                Ff(i,j) = (double)(F(i,j)+1);
            }
        }

        matrix_to_list(Vf,Vf_);
        matrix_to_list(Ff,Ff_);

        mxSetCell(hulls_V_ptr, num_levels-1-cur_num_levels, writeMatrix(Vf_));
        mxSetCell(hulls_F_ptr, num_levels-1-cur_num_levels, writeMatrix(Ff_));

        cur_num_levels++;

      }

    // Output to Matlab
//    plhs[0] = writeMatrix(Vf_);
//    plhs[1] = writeMatrix(Ff_);
    plhs[0] = hulls_V_ptr;
    plhs[1] = hulls_F_ptr;

    // delete stuff
    Vf_.clear();
    Ff_.clear();
    V0_.clear();
    F0_.clear();
    levels_.clear();

    // Restore the std stream buffer Important!
    std::cout.rdbuf(outbuf);
    return;

    }
