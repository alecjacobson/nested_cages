// Our header files
#include "io.h"
#include "cgal.h"
#include "flow.h"

// libigl includes
#include <igl/doublearea.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>

// useful namespaces
using namespace Eigen;
using namespace igl;

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
  using namespace igl::copyleft::cgal;

  if (argc==1){
    cout<<"Usage: ./nested_cages [filename.(off)] q L(1) L(2) ... L(k)"<<endl;
    cout<<"  q is the quadrature order for the shrinking flow"<<endl;
    cout<<"  L(1) L(2) ... L(k) is the number of faces for each cage"<<endl;
    return 0;
  }
  
  // number of layers
  int k = argc-4;
  cout << "number of layers = " << k << endl;

  // quadrature order
  int quad_order = atoi(argv[2]);

  Surface_mesh M; 
  std::ifstream is(argv[1]) ; is >> M ;

  // // output input mesh as level 0
  char* filename; 
  char *suffix;
  if (asprintf(&filename, "%s%s", argv[argc-1], "_0.off")!=-1){
    std::ofstream os( filename ) ; os << M;
  } else {
    cout << "unable to allocate space for output file name"  << endl;
    return 0;
  }

  // Converto CGAL's M to (V,F) format
  MatrixXd V0;
  MatrixXi F0;
  polyhedron_to_mesh(M,V0,F0); 

  // for now output CGAL decimations
  int L[k];
  bool adaptive = true;
  Surface_mesh M_hat;
  for(int i = 0;i<k;i++){
    std::ifstream is_file(argv[i+3]);
    // check if argv is a valid file
    if (is_file){
      is_file >> M_hat;
    } else{
      // first check if last charcater of argv[i+2] is r. If it is, drop the 'r' and adaptive = false
      adaptive = remove_all_chars_and_count(argv[i+3], 'r')==0;
      // else check if argv[i+2] is a valid integer (throw an error if it is not)
      if (!legal_int(argv[i+2])){
        cout << "you have to pass integer values or valid input deimatations"  << endl;
        cout << "the invalid argument you have passed is " << argv[i+3] << endl;
        return 0;
      }
      L[i] = atoi(argv[i+3]);
      float ratio = (1.*L[i])/(1.*M.size_of_facets()); 
      M_hat = M;
      decimate_CGAL(&M_hat,ratio,adaptive);

      // Check if decimations self-intersect. If they do, throw an error and quit 
      // (replace by Meshfix in the future)
	  MatrixXd V_coarse;
      MatrixXi F_coarse;
      polyhedron_to_mesh(M_hat,V_coarse,F_coarse); 
      RemeshSelfIntersectionsParam params;
      params.detect_only = true;
      params.first_only = true;
      MatrixXd tempV;
      MatrixXi tempF;
      MatrixXi IF;
      VectorXi J;
      VectorXi IM;
      remesh_self_intersections(V_coarse,F_coarse,params,tempV,tempF,IF,J,IM);
      if (IF.rows()>0){
      	cout << i+1 << "-th input decimation self-intersects. Quitting...  " << endl;  
      	return 0;
      }


	  // calculate triangle areas for initial mesh (will be used
	  // to define the integral at _every_ step - the metric is fixed)
	  VectorXd area_0;
	  doublearea(V0,F0,area_0);
	  area_0 = 0.5*area_0;
      // Precompute matrix that convert gradients at quadrature points to gradients at mesh vertices
  	  SparseMatrix<double> A_qv;
      gradQ_to_gradV(V0, F0, area_0, quad_order, A_qv);
  	  // Flow M inside M_hat
      MatrixXd V;
      flow_fine_inside_coarse(V0,F0,V_coarse,F_coarse,quad_order,A_qv,V);

      // Reinflate

      // Ouput cage
      if ((asprintf(&suffix,"_%d.off",i+1)!=-1) && (asprintf(&filename, "%s%s", argv[argc-1], suffix)!=-1)){
        std::ofstream os( filename ) ; os << M_hat;
      } else {
        cout << "unable to allocate space for output file name"  << endl;
        return 0;
      }

      M = M_hat;
      // back toi adaptive decimation (standard)
      adaptive = true;
    }
  }
  
  free(filename);
  return 1;
}
