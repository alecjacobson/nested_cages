#include "flow.h"

#include <igl/signed_distance.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>

#include <CGAL/Cartesian.h>

#include <iostream>


// convert gradient at quadrature points to gradient ate mesh vertices 
// (follow grad_quadrature_to_vertices.m)

void gradQ_to_gradV(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  const Eigen::MatrixXd & area_0, 
  const int quad_order,
  Eigen::SparseMatrix<double> & A)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;

  typedef Triplet<double> Tripletd;

  if (quad_order==1)
  {
    A.resize(V0.rows(),F0.rows());
    A.reserve(3*F0.rows());

    vector<Triplet<double> > IJV;
    IJV.reserve(3*F0.rows());

    for (int k=0;k<F0.rows();k++)
    {
      IJV.push_back(Tripletd(F0(k,0),k,area_0(k)/3.));
      IJV.push_back(Tripletd(F0(k,1),k,area_0(k)/3));
      IJV.push_back(Tripletd(F0(k,2),k,area_0(k)/3));
    }
    A.setFromTriplets(IJV.begin(),IJV.end());

  } else if (quad_order==2)
  {

    A.resize(V0.rows(),3*F0.rows());

    vector<Triplet<double> > IJV;
    IJV.reserve(6*F0.rows());

    for (int k=0;k<F0.rows();k++)
    {
      IJV.push_back(Tripletd(F0(k,0),k,area_0(k)/6));
      IJV.push_back(Tripletd(F0(k,0),2*F0.rows()+k,area_0(k)/6));

      IJV.push_back(Tripletd(F0(k,1),F0.rows()+k,area_0(k)/6));
      IJV.push_back(Tripletd(F0(k,1),k,area_0(k)/6));

      IJV.push_back(Tripletd(F0(k,2),2*F0.rows()+k,area_0(k)/6));
      IJV.push_back(Tripletd(F0(k,2),F0.rows()+k,area_0(k)/6));
    }
    A.setFromTriplets(IJV.begin(),IJV.end());

  } else if (quad_order==3)
  {

    A.resize(V0.rows(),4*F0.rows());

    vector<Triplet<double> > IJV;
    IJV.reserve(12*F0.rows());

    for (int k=0;k<F0.rows();k++)
    {
      IJV.push_back(Tripletd(F0(k,0),k,-(3.0/16.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,0),F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,0),2*F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,0),3*F0.rows()+k,(55.0/144.0)*area_0(k)));

      IJV.push_back(Tripletd(F0(k,1),k,-(3.0/16.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,1),F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,1),2*F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,1),3*F0.rows()+k,(55.0/144.0)*area_0(k)));

      IJV.push_back(Tripletd(F0(k,2),k,-(3.0/16.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,2),F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,2),2*F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,2),3*F0.rows()+k,(55.0/144.0)*area_0(k)));
    }
    A.setFromTriplets(IJV.begin(),IJV.end());

    } else 
    {
      assert(false && "quad_order should be 1, 2, or 3");
      cout << "gradQ_to_gradV, quadrature order " << quad_order << " is not implemented"  << endl;
    }


}

void signed_distance_direction(
  const Eigen::MatrixXd & P, 
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  Eigen::MatrixXd & D)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  MatrixXd C,N;
  VectorXi I;
  VectorXd S;
  igl::signed_distance(P,V,F,igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,S,I,C,N);
  // next: continue writing as in signed_distance_direction.m, case 3 
}


void grad_energy(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  int quad_order, 
  const Eigen::SparseMatrix<double> & A_qv, 
  Eigen::VectorXd & grad)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;

  if (quad_order==1)
  {
    cout << "grad_energy, quadrature order = " << quad_order  << endl;
    // quadrature points: barycenters of triangles
    MatrixXd p123(F.rows(),3);
    for (int k=0; k<F.rows(); k++)
    {
      p123(k,0) = (1.0/3.0)*(V(F(k,0),0)+V(F(k,1),0)+V(F(k,2),0));
      p123(k,1) = (1.0/3.0)*(V(F(k,0),1)+V(F(k,1),1)+V(F(k,2),1));
      p123(k,2) = (1.0/3.0)*(V(F(k,0),2)+V(F(k,1),2)+V(F(k,2),2));
    }
    MatrixXd quad_points = p123;

    MatrixXd grad_Q;
    signed_distance_direction(quad_points,V_coarse,F_coarse,grad_Q);

    // [S,I,C,N] = signed_distance(P,V,F,'SignedDistanceType','pseudonormal');

  } else if (quad_order==2)
  {
    cout << "grad_energy, quadrature order = " << quad_order  << endl;
  } else if (quad_order==3)
  {
    cout << "grad_energy, quadrature order = " << quad_order  << endl;
  } else
  {
    assert(false && "quad_order should be 1, 2, or 3");
    cout << "grad_energy, quadrature order " << quad_order << " is not implemented"  << endl;
  }

  grad.resize(V.rows());
  for (int k=0;k<grad.rows();k++)
  {
    grad(k) = 0.0;
  }

}

// one step function
void flow_one_step(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const int quad_order, 
  const Eigen::SparseMatrix<double> & A_qv, 
  Eigen::MatrixXd & V_new)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  VectorXd grad;
  grad_energy(V, F, V_coarse, F_coarse, quad_order, A_qv, grad);

  V_new = MatrixXd::Zero(V.rows(),3);
}


// shrink until inside function
void flow_fine_inside_coarse(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  int quad_order, 
  const Eigen::SparseMatrix<double> & A_qv,
  Eigen::MatrixXd & V)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;

  // while there are intersections or some negative winding number, keep flowing

  // **Alec: notice that we cannot pass V as input and output, instead wrap the 
  // input inside MatrixXd to force compiler to make a copy in this case.
  flow_one_step(MatrixXd(V), F0, V_coarse, F_coarse, quad_order, A_qv, V);
}
