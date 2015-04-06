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
#include <igl/read_triangle_mesh.h>
#include <igl/edges.h>
#include <igl/unique_edge_map.h>
#include <igl/unique.h>
#include <igl/intersect.h>
#include <igl/list_to_matrix.h>
#include <igl/Viewer/viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

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

  VectorXi EMAP;
  MatrixXi E,EF,EI;
  edge_flaps(F,E,EMAP,EF,EI);

  typedef std::set<std::pair<double,int> > PriorityQueue;
  PriorityQueue Q;
  std::vector<PriorityQueue::iterator > Qit(E.rows());
  for(int e = 0;e<E.rows();e++)
  {
    Qit[e] = Q.emplace((V.row(E(e,0))-V.row(E(e,1))).norm(),e).first;
  }

  // number of edges
  int num_collapsed = 0;

  const auto &pre_draw = [&](igl::Viewer & viewer)->bool
  {
    if(viewer.core.is_animating && !Q.empty())
    {
      // collapse edge 
      for(int j = 0;j<100;j++)
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
          collapse_edge(p.second,V,F,E,EMAP,EF,EI,e1,e2,f1,f2);
        if(collapsed)
        {
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
                Q.erase(Qit[e]);
                Qit[e] = 
                  Q.emplace((V.row(E(e,0))-V.row(E(e,1))).norm(),e).first;
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
