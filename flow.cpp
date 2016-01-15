#include "flow.h"

#include <igl/signed_distance.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/normalize_row_lengths.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/copyleft/cgal/intersect_other.h>
#include <igl/winding_number.h>

#include <CGAL/Cartesian.h>

#include <iostream>

// compute the matrix that converts gradient at quadrature 
// points to gradient at mesh vertices

void gradQ_to_gradV(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  const Eigen::VectorXd & area_0, 
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
      assert(F0(k,0)<V0.rows() && F0(k,1)<V0.rows() && F0(k,2)<V0.rows());
      assert(k<F0.rows());
      IJV.push_back(Tripletd(F0(k,0),k,area_0(k)/3));
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
      IJV.push_back(Tripletd(F0(k,1),F0.rows()+k,(55.0/144.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,1),2*F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,1),3*F0.rows()+k,(5.0/72.0)*area_0(k)));

      IJV.push_back(Tripletd(F0(k,2),k,-(3.0/16.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,2),F0.rows()+k,(5.0/72.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,2),2*F0.rows()+k,(55.0/144.0)*area_0(k)));
      IJV.push_back(Tripletd(F0(k,2),3*F0.rows()+k,(5.0/72.0)*area_0(k)));
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
  MatrixXd dif = C-P;
  igl::normalize_row_lengths(dif,D);
  for (int k=0; k<S.rows(); k++){
    if (S(k)<=-1e-5){
      D.row(k) = -D.row(k);
    }
    else if (S(k)<1e-5){
      D.row(k) = -N.row(k);
    }
  }
}


void grad_energy(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv,
  const Eigen::SparseMatrix<double> & M_inv, 
  Eigen::MatrixXd & grad)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;

  // Compute quad_order from A_qv
  int quad_order = [&F,&A_qv]()->int
  {
    if(A_qv.cols() == F.rows())
    {
      return 1;
    }else if(A_qv.cols() == 3*F.rows())
    {
      return 2;
    }else if(A_qv.cols() == 4*F.rows())
    {
      return 3;
    }else
    {
      assert(false && "#A_qv should be multiple of #F");
      return -1;
    }
  }();

  if (quad_order==1)
  {
    // quadrature points: barycenters of triangles
    MatrixXd p123(F.rows(),3);
    for (int k=0; k<F.rows(); k++)
    {
      p123(k,0) = (1.0/3.0)*(V(F(k,0),0)+V(F(k,1),0)+V(F(k,2),0));
      p123(k,1) = (1.0/3.0)*(V(F(k,0),1)+V(F(k,1),1)+V(F(k,2),1));
      p123(k,2) = (1.0/3.0)*(V(F(k,0),2)+V(F(k,1),2)+V(F(k,2),2));
    }
    MatrixXd quad_points = p123;

    MatrixXd grad_Q(F.rows(),3);
    signed_distance_direction(quad_points,V_coarse,F_coarse,grad_Q);
    grad_Q = -grad_Q;
    grad = M_inv*(A_qv*grad_Q);

  } else if (quad_order==2)
  {

    MatrixXd p12(F.rows(),3);
    MatrixXd p23(F.rows(),3);
    MatrixXd p31(F.rows(),3);
    for (int k=0; k<F.rows(); k++)
    {
      p12(k,0) = (1.0/2.0)*(V(F(k,0),0)+V(F(k,1),0));
      p12(k,1) = (1.0/2.0)*(V(F(k,0),1)+V(F(k,1),1));
      p12(k,2) = (1.0/2.0)*(V(F(k,0),2)+V(F(k,1),2));

      p23(k,0) = (1.0/2.0)*(V(F(k,1),0)+V(F(k,2),0));
      p23(k,1) = (1.0/2.0)*(V(F(k,1),1)+V(F(k,2),1));
      p23(k,2) = (1.0/2.0)*(V(F(k,1),2)+V(F(k,2),2));

      p31(k,0) = (1.0/2.0)*(V(F(k,2),0)+V(F(k,0),0));
      p31(k,1) = (1.0/2.0)*(V(F(k,2),1)+V(F(k,0),1));
      p31(k,2) = (1.0/2.0)*(V(F(k,2),2)+V(F(k,0),2));
    }

    // concatenate quadrature points
    MatrixXd quad_points(p12.rows()+p23.rows()+p31.rows(), p12.cols());
    quad_points << p12,
         p23,
         p31;

    MatrixXd grad_Q(3*F.rows(),3);
    signed_distance_direction(quad_points,V_coarse,F_coarse,grad_Q);
    grad_Q = -grad_Q;
    grad = M_inv*(A_qv*grad_Q);
    // the correct would be grad = M\(A_qv*grad_Q). Ask Alec waht solver he uses

  } else if (quad_order==3)
  {

    MatrixXd p1(F.rows(),3);
    MatrixXd p2(F.rows(),3);
    MatrixXd p3(F.rows(),3);
    MatrixXd p4(F.rows(),3);
    for (int k=0; k<F.rows(); k++)
    {
      p1(k,0) = (1.0/3.0)*(V(F(k,0),0)+V(F(k,1),0)+V(F(k,2),0));
      p1(k,1) = (1.0/3.0)*(V(F(k,0),1)+V(F(k,1),1)+V(F(k,2),1));
      p1(k,2) = (1.0/3.0)*(V(F(k,0),2)+V(F(k,1),2)+V(F(k,2),2));

      p2(k,0) = (2.0/15.0)*V(F(k,0),0) + (11.0/15.0)*V(F(k,1),0)+ (2.0/15.0)*V(F(k,2),0);
      p2(k,1) = (2.0/15.0)*V(F(k,0),1) + (11.0/15.0)*V(F(k,1),1)+ (2.0/15.0)*V(F(k,2),1);
      p2(k,2) = (2.0/15.0)*V(F(k,0),2) + (11.0/15.0)*V(F(k,1),2)+ (2.0/15.0)*V(F(k,2),2);

      p3(k,0) = (2.0/15.0)*V(F(k,0),0) + (2.0/15.0)*V(F(k,1),0)+ (11.0/15.0)*V(F(k,2),0);
      p3(k,1) = (2.0/15.0)*V(F(k,0),1) + (2.0/15.0)*V(F(k,1),1)+ (11.0/15.0)*V(F(k,2),1);
      p3(k,2) = (2.0/15.0)*V(F(k,0),2) + (2.0/15.0)*V(F(k,1),2)+ (11.0/15.0)*V(F(k,2),2);

      p4(k,0) = (11.0/15.0)*V(F(k,0),0) + (2.0/15.0)*V(F(k,1),0)+ (2.0/15.0)*V(F(k,2),0);
      p4(k,1) = (11.0/15.0)*V(F(k,0),1) + (2.0/15.0)*V(F(k,1),1)+ (2.0/15.0)*V(F(k,2),1);
      p4(k,2) = (11.0/15.0)*V(F(k,0),2) + (2.0/15.0)*V(F(k,1),2)+ (2.0/15.0)*V(F(k,2),2);
    }

    // concatenate quadrature points
    MatrixXd quad_points(p1.rows()+p2.rows()+p3.rows()+p4.rows(), p1.cols());
    quad_points << p1,
         p2,
         p3,
         p4;

    MatrixXd grad_Q(4*F.rows(),3);
    signed_distance_direction(quad_points,V_coarse,F_coarse,grad_Q);
    grad_Q = -grad_Q;
    grad = M_inv*(A_qv*grad_Q);
    // the correct would be grad = M\(A_qv*grad_Q). Ask Alec waht solver he uses

  } else
  {
    assert(false && "quad_order should be 1, 2, or 3");
    cout << "grad_energy, quadrature order " << quad_order << " is not implemented"  << endl;
  }

}

// one step function
void flow_one_step(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv, 
  const Eigen::SparseMatrix<double> & M,
  const double delta_t,  
  Eigen::MatrixXd & V_new)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  MatrixXd grad;
  grad_energy(V, F, V_coarse, F_coarse, A_qv, M, grad);

  V_new = V - delta_t*grad;
}

double diameter(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXd & V_coarse)
{

  using namespace Eigen;
  using namespace std;

  double min_x, min_y, min_z, max_x, max_y, max_z;

  VectorXd min_V0 = V0.colwise().minCoeff();
  VectorXd min_V_coarse = V_coarse.colwise().minCoeff();
  min_x = fmin(min_V0(0),min_V_coarse(0));
  min_y = fmin(min_V0(1),min_V_coarse(1));
  min_z = fmin(min_V0(2),min_V_coarse(2));

  VectorXd max_V0 = V0.colwise().maxCoeff();
  VectorXd max_V_coarse = V_coarse.colwise().maxCoeff();
  max_x = fmax(max_V0(0),max_V_coarse(0));
  max_y = fmax(max_V0(1),max_V_coarse(1));
  max_z = fmax(max_V0(2),max_V_coarse(2));

  float diam = fmax(max_x-min_x,max_y-min_y);
  diam = fmax(diam,max_z-min_z);

  return diam;
}

bool flow_fine_inside_coarse(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv,
  std::stack<Eigen::MatrixXd> & H)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  using namespace igl::copyleft::cgal;

  // while there are intersections or some negative winding number, keep flowing
  SparseMatrix<double> M;
  massmatrix(V0,F0,MASSMATRIX_TYPE_BARYCENTRIC,M);
  SparseMatrix<double> M_inv;
  invert_diag(M,M_inv);

  MatrixXd V = V0;
  H.push(V);
  // calculate diameter of the meshes to scale step size
  double diam = diameter(V0,V_coarse);
  double delta_t = diam*1e-3;
  MatrixXi IF;
  intersect_other(V,F0,V_coarse,F_coarse,true,IF);

  VectorXd W(1); // winding number of the first point
  winding_number(V_coarse,F_coarse,V.row(0),W);

  int step = 1; 

  while (IF.rows()>0 || W[0]<1e-10){
    #ifdef VERBOSE_DEBUG
      cout << "Flow step " << step << ":Coarse and fine mesh intersect " << endl;
    #endif
    flow_one_step(MatrixXd(V), F0, V_coarse, F_coarse, A_qv, M_inv, delta_t, V);
    intersect_other(V,F0,V_coarse,F_coarse,true,IF);
    step = step+1;
    // **Alec: notice that we cannot pass V as input and output, instead wrap the 
    // input inside MatrixXd to force compiler to make a copy in this case.
    winding_number(V_coarse,F_coarse,V.row(0),W);
    H.push(V);
    // Quit after 1000 steps of the flow and return false
    if (step>1000) return false; 
  }
  return true;
}
