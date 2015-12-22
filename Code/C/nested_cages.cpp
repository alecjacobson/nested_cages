// Our header files
#include "io.h"
#include "cgal.h"
#include "flow.h"
#include "reinflate.h"

// libigl includes
#include <igl/doublearea.h>
#include <igl/writeOFF.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <igl/copyleft/cgal/intersect_other.h>
 
// useful namespaces
using namespace Eigen;
using namespace igl;
using namespace std;

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
    cout << "Usage: ./nested_cages input.off q L(1) L(2) ... L(k) EnergyExpansion EnergyFinal output.off" << endl;
    cout << "" << endl;
    cout << "  q is the quadrature order for the shrinking flow" << endl;
    cout << "" << endl;
    cout << "  L(1) > L(2) > ... > L(k) is the number of faces for each cage" << endl;
    cout << "  If L(k) is followed by 'r' the initial decimation for this cage will be regular (adaptive if no 'r')" << endl;
    cout << "  Each L(k) can be replace by a filed with an input decimation " << endl;
    cout << "" << endl;
    cout << "  EnergyExpansion is the energy to be minimized for the re-inflation" << endl;
    cout << "  Energies implemented: None, DispStep, DispInitial, Volume, SurfARAP, VolARAP " << endl;
        cout << "" << endl;
    cout << "  EnergyFinal is the energy to be minimized after the re-inflation (additional processing)" << endl;
    cout << "  Energies implemented: None, DispStep, DispInitial, Volume, SurfARAP, VolARAP " << endl;
    return 0;
  }
  
  // number of layers
  int k = argc-6;
  cout << "number of layers = " << k << endl;

  // quadrature order
  int quad_order = atoi(argv[2]);

  // read input mesh
  Surface_mesh M; 
  std::ifstream is(argv[1]) ; is >> M ;

  // output input mesh as level output_0.off
  char* filename; 
  char *suffix;
  if (asprintf(&filename, "%s%s", argv[argc-1], "_0.off")!=-1){
    std::ofstream os( filename ) ; os << M;
  } else {
    cout << "unable to allocate space for output file name"  << endl;
    return 0;
  }

  // Converto CGAL's M to LibiIGL/Eigen (V,F) format
  MatrixXd V0;
  MatrixXi F0;
  polyhedron_to_mesh(M,V0,F0); 

  // First fine mesh is the input mesh
  MatrixXd V = V0;
  MatrixXi F = F0;

  // vector where each entry is the number of faces for eachj level
  int L[k];
  // standard CGAL decimation is adaptive
  bool adaptive = true;
  // declare input decimation
  Surface_mesh M_hat;

  // Loop over levels
  for(int i = 0;i<k;i++){
    // check if argv[i+3] is a valid file. If so, it will be used as input decimation 
    std::ifstream is_file(argv[i+3]);
    if (is_file)
    {
      is_file >> M_hat;
    } 
    // otherwise the input decimation is computed by decimating the previous layer
    // with CGAL decimation (regular or adaptive)
    else
    {
      // first check if last charcater of argv[i+2] is r. If it is, drop the 'r' and adaptive = false
      adaptive = remove_all_chars_and_count(argv[i+3], 'r')==0;
      // check if argv[i+2] is a valid integer (throw an error if it is not)
      if (!legal_int(argv[i+2])){
        cout << "you have to pass integer values or valid input deimatations"  << endl;
        cout << "the invalid argument you have passed is " << argv[i+3] << endl;
        return 0;
      }
      // specified number of faces for this cage
      L[i] = atoi(argv[i+3]);
      // value to pass to the decimator
      float ratio = (1.*L[i])/(1.*M.size_of_facets()); 
      // the previously computed will be decimated
      M_hat = M;
      // decimate
      decimate_CGAL(&M_hat,ratio,adaptive);
    }

    // Check if decimations self-intersect. If they do, throw an error and quit 
    // (replace by Meshfix in the future)
	  MatrixXd V_coarse;
    MatrixXi F_coarse;
    // Cnvert decimation to LibIGL/Eigen format
    polyhedron_to_mesh(M_hat,V_coarse,F_coarse); 
    // Parameters to call function to check for decimation's self-intersections
    RemeshSelfIntersectionsParam params;
    params.detect_only = true;
    params.first_only = true;
    MatrixXd tempV;
    MatrixXi tempF;
    MatrixXi IF;
    VectorXi J;
    VectorXi IM;
    remesh_self_intersections(V_coarse,F_coarse,params,tempV,tempF,IF,J,IM);
    // If input coarse mesh self-intersect, remove self-intersections with Meshfix (to-do)
    if (IF.rows()>0){
    	cout << i+1 << "-th input decimation self-intersects. Quitting...  " << endl;  
    	return 0;
    }

	  // calculate triangle areas for initial mesh (will be used
	  // to define the integral at _every_ step - the metric is fixed)
	  VectorXd area_0;
	  doublearea(V,F,area_0);
	  area_0 = 0.5*area_0;
    // Precompute matrix that convert gradients at quadrature points to gradients at mesh vertices
  	SparseMatrix<double> A_qv;
    gradQ_to_gradV(V, F, area_0, quad_order, A_qv);
  	// Flow M inside M_hat and save the result to a stack M of flow meshes
    stack<MatrixXd> H;
    flow_fine_inside_coarse(V,F,V_coarse,F_coarse,A_qv,H);

    // Reinflate and output to cage to C
    MatrixXd C;
    reinflate(H,F,V_coarse,F_coarse,argv[argc-3],argv[argc-2],C);

    // sanity check: cage should never self-intersect at this stage
    remesh_self_intersections(C,F_coarse,params,tempV,tempF,IF,J,IM);
    if (IF.rows()>0){
      cout << i+1 << "-th output cage self-intersects. ERROR! Quitting...  " << endl;  
      return 0;
    }

    // sanity check: cage should never intersect input decimation
    intersect_other(C,F_coarse,V,F,true,IF);
    if (IF.rows()>0){
      cout << i+1 << "-th output cage intersect previous cage. ERROR! Quitting...  " << endl;  
      return 0;
    }

    // output cage is the input for the next level
    M.clear();
    mesh_to_polyhedron(C,F_coarse,M); 
    V = C;
    F = F_coarse;

    // Output cage to file output_i.off
    if ((asprintf(&suffix,"_%d.off",i+1)!=-1) && (asprintf(&filename, "%s%s", argv[argc-1], suffix)!=-1)){
      writeOFF(filename,C,F_coarse);
    } else {
      cout << "unable to allocate space for output file name"  << endl;
      return 0;
    }
    // back to adaptive decimation (standard)
    adaptive = true;
  }
  
  free(filename);
  return 1;
}
