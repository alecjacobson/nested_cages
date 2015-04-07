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

#include "collapse_edge.h"
#include "edge_flaps.h"
//#include "linprog.h"
//#include "mosek_linprog.h"
#include "eiquadprog.hpp"
#include <igl/read_triangle_mesh.h>
#include <igl/centroid.h>
#include <igl/matlab_format.h>
#include <igl/edges.h>
#include <igl/unique_edge_map.h>
#include <igl/unique.h>
#include <igl/intersect.h>
#include <igl/list_to_matrix.h>
#include <igl/Viewer/viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>
#include <algorithm>

// For debuggin'
int at(
  Eigen::MatrixXi & M,
  const int i,
  const int j)
{
  return M(i,j);
}

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  string
    filename("/Users/ajx/Dropbox/models/basic shapes/sphere-hires.obj");
    //filename("/usr/local/igl/libigl/examples/shared/cheburashka.off");
    //filename("/usr/local/igl/libigl/examples/shared/decimated-knight.obj");
  if(argc>=2)
  {
    filename = argv[1];
  }
  MatrixXd V,OV;
  MatrixXi F,OF;
  read_triangle_mesh(filename,OV,OF);
  igl::Viewer viewer;

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

  //// variables for mosek task, env and result code
  //MSKenv_t env;
  //// Create the MOSEK environment
  //mosek_guarded(MSK_makeenv(&env,NULL));
  //// initialize mosek environment
  //mosek_guarded(MSK_initenv(env));

  const auto & minimal_volume_and_outer_hull = [
    //&env
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
      const double w = 1;
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
      Qit[e] = Q.emplace(cost,e).first;
    }
    num_collapsed = 0;
    viewer.data.clear();
    viewer.data.set_mesh(V,F);
    viewer.data.set_face_based(true);
  };

  const auto &pre_draw = [&](igl::Viewer & viewer)->bool
  {
    if(viewer.core.is_animating && !Q.empty())
    {
      bool something_collapsed = false;
      // collapse edge 
      const int max_iter = std::ceil(0.01*Q.size());
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
            if(F(n,0) != COLLAPSE_EDGE_NULL ||
               F(n,1) != COLLAPSE_EDGE_NULL ||
               F(n,2) != COLLAPSE_EDGE_NULL)
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
                Qit[e] = Q.emplace(cost,e).first;
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
      if(something_collapsed)
      {
        viewer.data.clear();
        viewer.data.set_mesh(V,F);
        viewer.data.set_face_based(true);
      }
    }
    return false;
  };

  const auto &key_down = [&](igl::Viewer &viewer,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case ' ':
        viewer.core.is_animating ^= 1;
        break;
      case 'R':
      case 'r':
        reset();
        break;
      default:
        return false;
    }
    return true;
  };

  reset();
  viewer.core.is_animating = true;
  viewer.callback_key_down = key_down;
  viewer.callback_pre_draw = pre_draw;
  //MSK_deleteenv(&env)
  return viewer.launch();
}
