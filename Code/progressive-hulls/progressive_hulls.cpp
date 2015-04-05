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
#rm -f $TEMP
# END BASH SCRIPT
exit
*/

#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/edges.h>
#include <igl/unique_edge_map.h>
#include <igl/unique.h>
#include <igl/intersect.h>
#include <igl/list_to_matrix.h>
#include <igl/Viewer/viewer.h>
#include <iostream>

int at(
  Eigen::MatrixXi & M,
  const int i,
  const int j)
{
  return M(i,j);
}

// Assumes (V,F) is a closed manifold mesh (except for previouslly collapsed
// faces).
//
// Inputs:
//   e  index into E of edge to try to collapse
// Outputs:
//
// Collapses exactly two faces and exactly 3 edges from E (e and one side of
// each face gets collapsed to the other).
//
// TODO: Shouldn't need V2F or need to maintain V2F
// returns true if edge was collapsed
bool collapse_edge(
  const int e,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI,
  std::vector<std::vector<int> > & V2F)
{
  // Assign this to 0 rather than, say, -1 so that deleted elements will get
  // draw as degenerate elements at vertex 0 (which should always exist and
  // never get collapsed to anything else since it is the smallest index)
#define COLLAPSE_EDGE_NULL 0
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  const int eflip = E(e,0)>E(e,1);
  // source and destination
  const int s = eflip?E(e,1):E(e,0);
  const int d = eflip?E(e,0):E(e,1);
  if(s == COLLAPSE_EDGE_NULL && d==COLLAPSE_EDGE_NULL)
  {
    return false;
  }
  // check if edge collapse is valid: intersection of vertex neighbors of s and
  // d should be exactly 2+(s,d) = 4
  // http://stackoverflow.com/a/27049418/148668
  {
    // all neighbors (including self)
    const auto neighbors = [](
        const int v,
        const Eigen::MatrixXi & F,
        const vector<vector<int> > & V2F)->
        VectorXi
    {
      vector<int> N,uN;
      for(auto f : V2F[v])
      {
        N.push_back(F(f,0));
        N.push_back(F(f,1));
        N.push_back(F(f,2));
      }
      vector<size_t> _1,_2;
      igl::unique(N,uN,_1,_2);
      VectorXi uNm;
      list_to_matrix(uN,uNm);
      return uNm;
    };
    VectorXi Ns = neighbors(s,F,V2F);
    VectorXi Nd = neighbors(d,F,V2F);
    VectorXi N = igl::intersect(Ns,Nd);
    if(N.size() != 4)
    {
      return false;
    }
  }


  assert(s<d);
  // move source and destination to midpoint
  const RowVectorXd mid = (V.row(s) + V.row(d))*0.5;
  V.row(s) = mid;
  V.row(d) = mid;
  // merge face lists except edge's flaps
  {
    vector<int> faces;
    faces.reserve( V2F[s].size() + V2F[d].size() ); // preallocate memory
    for(auto p : {s,d})
    {
      for(auto f : V2F[p])
      {
        if(f != EF(e,0) && f != EF(e,1))
        {
          faces.push_back(f);
        }
      }
    }
    vector<size_t> _1,_2;
    igl::unique(faces,V2F[s],_1,_2);
  }

  const auto & kill_edge = [&E,&EI,&EF](const int e)
  {
    E(e,0) = COLLAPSE_EDGE_NULL;
    E(e,1) = COLLAPSE_EDGE_NULL;
    EF(e,0) = COLLAPSE_EDGE_NULL;
    EF(e,1) = COLLAPSE_EDGE_NULL;
    EI(e,0) = COLLAPSE_EDGE_NULL;
    EI(e,1) = COLLAPSE_EDGE_NULL;
  };

  // update edge info
  //cout<<"E(e,:): "<<E(e,0)<<" "<<E(e,1)<<endl;
  //cout<<"s: "<<s<<endl;
  //cout<<"d: "<<d<<endl;
  //cout<<"eflip: "<<eflip<<endl;
  //cout<<"F("<<EF(e,0)<<",:): "<<F(EF(e,0),0)<<" "<<F(EF(e,0),1)<<" "<<F(EF(e,0),2)<<endl;
  //cout<<"F("<<EF(e,1)<<",:): "<<F(EF(e,1),0)<<" "<<F(EF(e,1),1)<<" "<<F(EF(e,1),2)<<endl;
  // for each flap
  const int m = F.rows();
  for(int side = 0;side<2;side++)
  {
    const int f = EF(e,side);
    const int v = EI(e,side);
    const int sign = (eflip==0?1:-1)*(1-2*side);
    // next edge emanating from d
    //cout<<"  sign: "<<sign<<endl;
    //cout<<"  f: "<<f<<endl;
    //cout<<"  v: "<<v<<endl;
    //cout<<"  F(f,:): "<<F(f,0)<<" "<<F(f,1)<<" "<<F(f,2)<<endl;
    const int e1 = EMAP(f+m*((v+sign*1+3)%3));
    //cout<<"  e_1: "<<E(e1,0)<<" "<<E(e1,1)<<endl;
    // prev edge pointing to s
    const int e2 = EMAP(f+m*((v+sign*2+3)%3));
    //cout<<"  e_2: "<<E(e2,0)<<" "<<E(e2,1)<<endl;
    assert(E(e1,0) == d || E(e1,1) == d);
    assert(E(e2,0) == s || E(e2,1) == s);
    // face adjacent to f on e1, also incident on d
    const bool flip1 = EF(e1,1)==f;
    const int f1 = flip1 ? EF(e1,0) : EF(e1,1);
    //cout<<"  f1: "<<f1<<endl;
    assert(f1!=f);
    assert(F(f1,0)==d || F(f1,1)==d || F(f1,2) == d);
    // across from which vertex of f1 does e1 appear?
    const int v1 = flip1 ? EI(e1,0) : EI(e1,1);
    // Kill e1
    kill_edge(e1);
    // Kill f
    F(f,0) = COLLAPSE_EDGE_NULL;
    F(f,1) = COLLAPSE_EDGE_NULL;
    F(f,2) = COLLAPSE_EDGE_NULL;
    // map f1's edge on e1 to e2
    assert(EMAP(f1+m*v1) == e1);
    EMAP(f1+m*v1) = e2;
    // side opposite f2, the face adjacent to f on e2, also incident on s
    const int opp2 = (EF(e2,0)==f?0:1);
    assert(EF(e2,opp2) == f);
    EF(e2,opp2) = f1;
    EI(e2,opp2) = v1;
    // remap e2 from d to s
    E(e2,0) = E(e2,0)==d ? s : E(e2,0);
    E(e2,1) = E(e2,1)==d ? s : E(e2,1);
    //cout<<endl;
  }

  // finally, reindex faces and edges incident on d. Do this last so asserts
  // make sense.
  //cout<<"d: "<<d<<endl;
  for(auto f : V2F[d])
  {
    for(int v = 0;v<3;v++)
    {
      if(F(f,v) == d)
      {
        //cout<<"f: "<<f<<endl;
        //cout<<"v: "<<v<<endl;
        //cout<<"F("<<f<<",:): "<<F.row(f)<<endl;
        //cout<<"E("<<EMAP(f+m*((v+1)%3))<<",:)"<<E.row(EMAP(f+m*((v+1)%3)))<<endl;
        //cout<<"E("<<EMAP(f+m*((v+2)%3))<<",:)"<<E.row(EMAP(f+m*((v+2)%3)))<<endl;
        const int flip1 = (EF(EMAP(f+m*((v+1)%3)),0)==f)?1:0;
        const int flip2 = (EF(EMAP(f+m*((v+2)%3)),0)==f)?0:1;
        //cout<<"flip1: "<<flip1<<endl;
        //cout<<"flip2: "<<flip2<<endl;
        assert(
          E(EMAP(f+m*((v+1)%3)),flip1) == d ||
          E(EMAP(f+m*((v+1)%3)),flip1) == s);
        E(EMAP(f+m*((v+1)%3)),flip1) = s;
        assert(
          E(EMAP(f+m*((v+2)%3)),flip2) == d ||
          E(EMAP(f+m*((v+2)%3)),flip2) == s);
        E(EMAP(f+m*((v+2)%3)),flip2) = s;
        F(f,v) = s;
        break;
      }
    }
  }
  //cout<<"-----------------------------"<<endl<<endl;
  // Finally, "remove" this edge and its information
  kill_edge(e);

  return true;
}

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  string
    filename("/usr/local/igl/libigl/examples/shared/cheburashka.off");
  if(argc>=2)
  {
    filename = argv[1];
  }
  MatrixXd V;
  MatrixXi F;
  read_triangle_mesh(filename,V,F);
  igl::Viewer viewer;
  viewer.data.set_mesh(V,F);
  viewer.data.set_face_based(true);
  MatrixXi allE,E;
  VectorXi EMAP;
  vector<vector<int> > uE2E;
  unique_edge_map(F,allE,E,EMAP,uE2E);
  //cout<<"F: "<<F.rows()<<" "<<F.cols()<<endl;
  //cout<<"EMAP: "<<EMAP.rows()<<" "<<EMAP.cols()<<endl;
  // EF  #EF by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
  //   F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
  //   e=(j->i)
  MatrixXi EF(E.rows(),2);
  MatrixXi EI(E.rows(),2);
  vector<vector<int> > V2F(V.rows());
  for(int f = 0;f<F.rows();f++)
  {
    for(int v = 0;v<3;v++)
    {
      // add to vertex-to-faces list
      V2F[F(f,v)].push_back(f);
      // get edge id
      const int e = EMAP(v*F.rows()+f);
      // See if this is left or right flap w.r.t. edge orientation
      if( F(f,(v+1)%3) == E(e,0) && F(f,(v+2)%3) == E(e,1))
      {
        EF(e,0) = f;
        EI(e,0) = v;
      }else
      {
        assert(F(f,(v+1)%3) == E(e,1) && F(f,(v+2)%3) == E(e,0));
        EF(e,1) = f;
        EI(e,1) = v;
      }
    }
  }


  // number of edges
  int num_collapsed = 0;
  int e = 0;

  const auto &pre_draw = [&](igl::Viewer & viewer)->bool
  {
    if(viewer.core.is_animating)
    {
      // collapse edge 
      for(int i = 0;i<100;i++)
      for(int j = 0;j<100;j++)
      {
        const bool collapsed =
          collapse_edge(e,V,F,E,EMAP,EF,EI,V2F);
        e = (e+1)%E.rows();
        if(collapsed)
        {
          num_collapsed++;
          break;
        }
      }
      viewer.data.clear();
      viewer.data.set_mesh(V,F);
      viewer.data.set_face_based(true);
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
      default:
        return false;
    }
    return true;
  };

  viewer.core.is_animating = true;
  viewer.callback_key_down = key_down;
  viewer.callback_pre_draw = pre_draw;
  return viewer.launch();
}
