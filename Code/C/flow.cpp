#include "flow.h"

#include <CGAL/Cartesian.h>


// convert gradient at quadrature points to gradient ate mesh vertices 
// (follow grad_quadrature_to_vertices.m)

SparseMatrix<double> gradQ_to_gradV(MatrixXd V0, MatrixXi F0, MatrixXd area_0, int quad_order){

	SparseMatrix<double> A;

	if (quad_order==1) {
		A.resize(V0.rows(),F0.rows());
		A.reserve(3*F0.rows());

		vector<Triplet<double> > IJV;
  		IJV.reserve(3*F0.rows());

  		for (int k=0;k<F0.rows();k++){
  			IJV.push_back(Triplet<double>(F0(k,0),k,area_0(k)/3));
  			IJV.push_back(Triplet<double>(F0(k,1),k,area_0(k)/3));
  			IJV.push_back(Triplet<double>(F0(k,2),k,area_0(k)/3));
  		}
  		A.setFromTriplets(IJV.begin(),IJV.end());

	} else if (quad_order==2){

		A.resize(V0.rows(),3*F0.rows());

		vector<Triplet<double> > IJV;
  		IJV.reserve(6*F0.rows());

  		for (int k=0;k<F0.rows();k++){
  			IJV.push_back(Triplet<double>(F0(k,0),k,area_0(k)/6));
  			IJV.push_back(Triplet<double>(F0(k,0),2*F0.rows()+k,area_0(k)/6));

  			IJV.push_back(Triplet<double>(F0(k,1),F0.rows()+k,area_0(k)/6));
  			IJV.push_back(Triplet<double>(F0(k,1),k,area_0(k)/6));

  			IJV.push_back(Triplet<double>(F0(k,2),2*F0.rows()+k,area_0(k)/6));
  			IJV.push_back(Triplet<double>(F0(k,2),F0.rows()+k,area_0(k)/6));
  		}
  		A.setFromTriplets(IJV.begin(),IJV.end());

	} else if (quad_order==3){

		A.resize(V0.rows(),4*F0.rows());

		vector<Triplet<double> > IJV;
  		IJV.reserve(12*F0.rows());

  		for (int k=0;k<F0.rows();k++){
  			IJV.push_back(Triplet<double>(F0(k,0),k,-(3.0/16.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,0),F0.rows()+k,(5.0/72.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,0),2*F0.rows()+k,(5.0/72.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,0),3*F0.rows()+k,(55.0/144.0)*area_0(k)));

  			IJV.push_back(Triplet<double>(F0(k,1),k,-(3.0/16.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,1),F0.rows()+k,(5.0/72.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,1),2*F0.rows()+k,(5.0/72.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,1),3*F0.rows()+k,(55.0/144.0)*area_0(k)));

  			IJV.push_back(Triplet<double>(F0(k,2),k,-(3.0/16.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,2),F0.rows()+k,(5.0/72.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,2),2*F0.rows()+k,(5.0/72.0)*area_0(k)));
  			IJV.push_back(Triplet<double>(F0(k,2),3*F0.rows()+k,(55.0/144.0)*area_0(k)));
  		}
  		A.setFromTriplets(IJV.begin(),IJV.end());

		} else {
			cout << "gradQ_to_gradV, quadrature order " << quad_order << " is not implemented"  << endl;
		}

	return A;

}

MatrixXd signed_distance_direction(MatrixXd P, MatrixXd V, MatrixXi F){
	
	MatrixXd C,N;
	VectorXi I;
  	VectorXd S;
  igl::signed_distance(P,V,F,igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,S,I,C,N);

  	// next: continue writing as in signed_distance_direction.m, case 3 


	MatrixXd D;

	


	return D;

}


// gradient of energy function
VectorXd grad_energy(MatrixXd V, MatrixXi F, MatrixXd V_coarse, MatrixXi F_coarse, int quad_order, SparseMatrix<double> A_qv){

	if (quad_order==1){
		cout << "grad_energy, quadrature order = " << quad_order  << endl;
		// quadrature points: barycenters of triangles
		MatrixXd p123(F.rows(),3);
		for (int k=0; k<F.rows(); k++){
			p123(k,0) = (1.0/3.0)*(V(F(k,0),0)+V(F(k,1),0)+V(F(k,2),0));
			p123(k,1) = (1.0/3.0)*(V(F(k,0),1)+V(F(k,1),1)+V(F(k,2),1));
			p123(k,2) = (1.0/3.0)*(V(F(k,0),2)+V(F(k,1),2)+V(F(k,2),2));
		}
		MatrixXd quad_points = p123;

		MatrixXd grad_Q;
  		grad_Q = signed_distance_direction(quad_points,V_coarse,F_coarse);

		// [S,I,C,N] = signed_distance(P,V,F,'SignedDistanceType','pseudonormal');

	} else if (quad_order==2){
		cout << "grad_energy, quadrature order = " << quad_order  << endl;
	} else if (quad_order==3){
		cout << "grad_energy, quadrature order = " << quad_order  << endl;
	} else{
		cout << "grad_energy, quadrature order " << quad_order << " is not implemented"  << endl;
	}

	VectorXd grad(V.rows());
	for (int k=0;k<grad.rows();k++){
		grad(k) = 0.0;
	}

	return grad; 

}

// one step function
MatrixXd flow_one_step(MatrixXd V, MatrixXi F, MatrixXd V_coarse, MatrixXi F_coarse, int quad_order, SparseMatrix<double> A_qv){


	VectorXd grad = grad_energy(V, F, V_coarse, F_coarse, quad_order, A_qv);

	MatrixXd V_new = MatrixXd::Zero(V.rows(),3);

	return V_new;

}

// shrink until inside function
void flow_fine_inside_coarse(MatrixXd V0, MatrixXi F0, MatrixXd V_coarse, MatrixXi F_coarse, int quad_order, SparseMatrix<double> A_qv){

  MatrixXd V = V0;
  MatrixXi F = F0;

  // while there are intersections or some negative winding number, keep flowing
  V = flow_one_step(V, F, V_coarse, F_coarse, quad_order, A_qv);

}
