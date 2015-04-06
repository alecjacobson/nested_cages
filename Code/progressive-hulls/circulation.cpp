#include "circulation.h"
#include <igl/list_to_matrix.h>

std::vector<int> circulation(
  const int e,
  const bool ccw,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI)
{
  // prepare output
  std::vector<int> N;
  N.reserve(6);
  const int m = F.rows();
  const auto & step = [&](
    const int e, 
    const int ff,
    int & ne, 
    int & nf)
  {
    assert((EF(e,1) == ff || EF(e,0) == ff) && "e should touch ff");
    const int fside = EF(e,1)==ff?1:0;
    const int nside = EF(e,0)==ff?1:0;
    const int nv = EI(e,nside);
    // get next face
    nf = EF(e,nside);
    // get next edge 
    const int dir = ccw?-1:1;
    ne = EMAP(nf+m*((nv+dir+3)%3));
  };
  // Always start with first face (ccw in step will be sure to turn right
  // direction)
  const int f0 = EF(e,0);
  int fi = f0;
  int ei = e;
  while(true)
  {
    N.push_back(fi);
    step(ei,fi,ei,fi);
    // back to start?
    if(fi == f0)
    {
      assert(ei == e);
      break;
    }
  }
  return N;
}

void circulation(
  const int e,
  const bool ccw,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI,
  Eigen::VectorXi & vN)
{
  std::vector<int> N = circulation(e,ccw,F,E,EMAP,EF,EI);
  igl::list_to_matrix(N,vN);
}
