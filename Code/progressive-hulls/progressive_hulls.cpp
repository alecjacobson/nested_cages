
#include "eiquadprog.hpp"
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
//#include <igl/linprog.h>
//#include <igl/mosek_linprog.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/boundary_facets.h>
#include <igl/circulation.h>
#include <igl/centroid.h>
#include <igl/matlab_format.h>
#include <igl/edges.h>
#include <igl/unique_edge_map.h>
#include <igl/unique.h>
#include <igl/writeOBJ.h>
#include <igl/intersect.h>
#include <igl/list_to_matrix.h>
#include <igl/Viewer/viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>
#include <algorithm>
#include <cstdlib>

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
  cout<<"Usage: ./progressive_hulls [filename.(off|obj|ply)] [w]"<<endl;
  cout<<"  [space]  toggle animation."<<endl;
  cout<<"  'r'  reset."<<endl;
  string filename("fandisk.off");
  double w = 0.5;
  switch(argc)
  {
    default:
    case 3:
      w = strtod(argv[2],NULL);
    case 2:
      filename = argv[1];
    case 1:
      break;
  }

  MatrixXd V,OV,U;
  MatrixXi F,OF,G;
  read_triangle_mesh(filename,OV,OF);
  igl::viewer::Viewer viewer;

  if(OF.size() == 0)
  {
    cerr<<"Input is empty."<<endl;
    return false;
  }
  vector<vector<vector<int > > > TT;
  vector<vector<vector<int > > > TTi;
  cout<<"triangle_triangle_adjacency..."<<endl;
  triangle_triangle_adjacency(OF,TT,TTi);

  VectorXi _1;
  if(!is_edge_manifold(OV,OF) || !is_vertex_manifold(OF,_1))
  {
    cerr<<"Input is not manifold."<<endl;
    return false;
  }
  if(boundary_facets<MatrixXi,MatrixXi>(OF).size() != 0)
  {
    cerr<<"Input is not closed."<<endl;
    return false;
  }

  VectorXi EMAP;
  MatrixXi E,EF,EI;
  typedef std::set<std::pair<double,int> > PriorityQueue;
  PriorityQueue Q;
  std::vector<PriorityQueue::iterator > Qit;
  // If an edge were collapsed, we'd collapse it to these points:
  MatrixXd C;
  int num_collapsed;
  int m;

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

  cost_and_placement = shortest_edge_and_midpoint;
  //cost_and_placement = minimal_volume_and_outer_hull; 

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
    m = F.rows();
    viewer.data.clear();
    viewer.data.set_mesh(V,F);
    viewer.data.set_face_based(true);
  };

  const auto &pre_draw = [&](igl::viewer::Viewer & viewer)->bool
  {
    if(viewer.core.is_animating && !Q.empty())
    {
      bool something_collapsed = false;
      // collapse edge 
      const int max_iter = std::ceil(0.01*Q.size());
      for(int j = 0;j<max_iter;j++)
      {
        if(!collapse_edge(cost_and_placement,V,F,E,EMAP,EF,EI,Q,Qit,C))
        {
          break;
        }
        something_collapsed = true;
        num_collapsed++;
        m-=2;
      }

      if(something_collapsed)
      {
        cout<<"m: "<<m<<endl;
        viewer.data.clear();
        viewer.data.set_mesh(V,F);
        viewer.data.set_face_based(true);
      }
    }
    return false;
  };

  const auto &key_down = [&](igl::viewer::Viewer &viewer,unsigned char key,int mod)->bool
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
