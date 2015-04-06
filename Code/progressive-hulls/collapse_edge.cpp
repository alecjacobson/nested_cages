#include "collapse_edge.h"
#include "circulation.h"
#include <igl/list_to_matrix.h>
#include <igl/unique.h>
#include <igl/intersect.h>
#include <vector>

bool collapse_edge(
  const int e,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI,
  int & a_e1,
  int & a_e2,
  int & a_f1,
  int & a_f2)
{
  // Assign this to 0 rather than, say, -1 so that deleted elements will get
  // draw as degenerate elements at vertex 0 (which should always exist and
  // never get collapsed to anything else since it is the smallest index)
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
    // all vertex neighbors around edge, including the two vertices of the edge
    const auto neighbors = [](
      const int e,
      const bool ccw,
      const Eigen::MatrixXi & F,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & EMAP,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI) 
    {
      vector<int> N,uN;
      vector<int> V2Fe = circulation(e, ccw,F,E,EMAP,EF,EI);
      for(auto f : V2Fe)
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
    VectorXi Ns = neighbors(e, eflip,F,E,EMAP,EF,EI);
    VectorXi Nd = neighbors(e,!eflip,F,E,EMAP,EF,EI);
    VectorXi N = igl::intersect(Ns,Nd);
    if(N.size() != 4)
    {
      return false;
    }
  }

  // Important to grab neighbors of d before monkeying with edges
  const std::vector<int> nV2Fd = circulation(e,!eflip,F,E,EMAP,EF,EI);

  // The following implementation strongly relies on s<d
  assert(s<d && "s should be less than d");
  // move source and destination to midpoint
  const RowVectorXd mid = (V.row(s) + V.row(d))*0.5;
  V.row(s) = mid;
  V.row(d) = mid;

  // Helper function to replace edge and associate information with NULL
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
  // for each flap
  const int m = F.rows();
  for(int side = 0;side<2;side++)
  {
    const int f = EF(e,side);
    const int v = EI(e,side);
    const int sign = (eflip==0?1:-1)*(1-2*side);
    // next edge emanating from d
    const int e1 = EMAP(f+m*((v+sign*1+3)%3));
    // prev edge pointing to s
    const int e2 = EMAP(f+m*((v+sign*2+3)%3));
    assert(E(e1,0) == d || E(e1,1) == d);
    assert(E(e2,0) == s || E(e2,1) == s);
    // face adjacent to f on e1, also incident on d
    const bool flip1 = EF(e1,1)==f;
    const int f1 = flip1 ? EF(e1,0) : EF(e1,1);
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
    if(side==0)
    {
      a_e1 = e1;
      a_f1 = f;
    }else
    {
      a_e2 = e1;
      a_f2 = f;
    }
  }

  // finally, reindex faces and edges incident on d. Do this last so asserts
  // make sense.
  //
  // Could actually skip first two, since those are always the two collpased
  // faces.
  for(auto f : nV2Fd)
  {
    for(int v = 0;v<3;v++)
    {
      if(F(f,v) == d)
      {
        const int flip1 = (EF(EMAP(f+m*((v+1)%3)),0)==f)?1:0;
        const int flip2 = (EF(EMAP(f+m*((v+2)%3)),0)==f)?0:1;
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
  // Finally, "remove" this edge and its information
  kill_edge(e);

  return true;
}

bool collapse_edge(
  const int e,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI)
{
  int e1,e2,f1,f2;
  return collapse_edge(e,V,F,E,EMAP,EF,EI,e1,e2,f1,f2);
}
