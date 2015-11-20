#include "flow.h"

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

void signed_distance(MatrixXd P, MatrixXd V, MatrixXi F, MatrixXd* C, MatrixXd* N, VectorXi* I, VectorXd* S){


	typedef CGAL::Simple_cartesian<double> Kernel;
	typedef 
	  CGAL::AABB_tree<
	  CGAL::AABB_traits<Kernel, 
	    CGAL::AABB_triangle_primitive<Kernel, 
	      std::vector<CGAL::Triangle_3<Kernel> >::iterator
	    > > > 
	  Tree;
	static MatrixXd g_V;
	static MatrixXi g_F;
	static SignedDistanceType g_sign_type = NUM_SIGNED_DISTANCE_TYPE;
	static Tree g_tree;
	static vector<CGAL::Triangle_3<Kernel> > g_T;
	static WindingNumberAABB<Eigen::Vector3d> g_hier;
	static MatrixXd g_FN,g_VN,g_EN;
	static MatrixXi g_E;
	static VectorXi g_EMAP;

	g_V = V;
    g_F = F;
    g_sign_type = SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
    g_tree.clear();
    g_T.clear();

    // Prepare distance computation
    point_mesh_squared_distance_precompute(V,F,g_tree,g_T);

    // "Signed Distance Computation Using the Angle Weighted Pseudonormal"
    // [Bærentzen & Aanæs 2005]
    per_face_normals(V,F,g_FN);
    per_vertex_normals(V,F,PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,g_FN,g_VN);
    per_edge_normals(V,F,PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,g_FN,g_EN,g_E,g_EMAP);

    (*N).resize(P.rows(),3);
  	(*S).resize(P.rows(),1);
  	(*I).resize(P.rows(),1);
  	(*C).resize(P.rows(),3);

  	for(int p = 0;p<P.rows();p++){
	    typedef typename Kernel::FT FT;
	    typedef typename Kernel::Point_3 Point_3;
	    typedef typename CGAL::Triangle_3<Kernel> Triangle_3;
	    typedef typename std::vector<Triangle_3>::iterator Iterator;
	    typedef typename CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
	    typedef typename CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
	    typedef typename CGAL::AABB_tree<AABB_triangle_traits> Tree;
	    typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;
	    const Point_3 q(P(p,0),P(p,1),P(p,2));
	    double s,sqrd;
	    Point_and_primitive_id pp;
	    Vector3d n(0,0,0);
	    signed_distance_pseudonormal(g_tree,g_T,g_F,g_FN,g_VN,g_EN,g_EMAP,q,s,sqrd,pp,n);
	    (*N).row(p) = n;
	    (*I)(p) = pp.second - g_T.begin();
	    (*S)(p) = s*sqrt(sqrd);
	    (*C)(p,0) = pp.first.x();
	    (*C)(p,1) = pp.first.y();
	    (*C)(p,2) = pp.first.z();
  	}

	return;

}

MatrixXd signed_distance_direction(MatrixXd P, MatrixXd V, MatrixXi F){
	
	MatrixXd C,N;
	VectorXi I;
  	VectorXd S;
  	signed_distance(P,V,F,&C,&N,&I,&S);

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