#include "nested_cages_lib.h"

// Meshfix include
#if WITH_MESHFIX
#include <meshfix.h>
#include <meshfix_eigen.h>
#endif

// Our header files
#include "cgal.h"
#include "flow.h"
#include "reinflate.h"

// libigl includes
#include <igl/doublearea.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <igl/copyleft/cgal/intersect_other.h>

#include <iostream>
#include <stack>

// mesh-in, mesh-out wrapper around MeshFix
static void meshfix_wrapper(
  const Eigen::MatrixXd & Vin,
  const Eigen::MatrixXi & Fin,
  Eigen::MatrixXd & Vout,
  Eigen::MatrixXi & Fout)
{
#if WITH_MESHFIX
  /////////////////////////////////////////////////////////////////////////
  // Convert to meshfix type, call meshfix, convert back from meshfix type
  T_MESH::TMesh::init(); // This is mandatory
  T_MESH::Basic_TMesh tin;
  meshfix_from_eigen_matrices(Vin,Fin,tin);
  meshfix(false,tin);
  meshfix_to_eigen_matrices(tin,Vout,Fout);
  /////////////////////////////////////////////////////////////////////////
#else
  (void)Vin; (void)Fin; (void)Vout; (void)Fout;
#endif
}

bool nested_cages(
  const Eigen::MatrixXd & V0,
  const Eigen::MatrixXi & F0,
  const int quad_order,
  const std::vector<int> & num_faces,
  const std::vector<bool> & regular,
  const std::vector<Eigen::MatrixXd> & V_decim,
  const std::vector<Eigen::MatrixXi> & F_decim,
  const std::string & energy_expansion,
  const std::string & energy_final,
  std::vector<Eigen::MatrixXd> & Vout,
  std::vector<Eigen::MatrixXi> & Fout,
  const std::function<void(int, const Eigen::MatrixXd &, const Eigen::MatrixXi &)>
    & per_level_callback)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  using namespace igl::copyleft::cgal;

  // number of layers
  const int k = (int)num_faces.size();

  Vout.clear();
  Fout.clear();
  Vout.reserve(k);
  Fout.reserve(k);

  // convert input mesh to CGAL format
  Surface_mesh M;
  mesh_to_polyhedron(V0,F0,M);

  // First fine mesh is the input mesh
  MatrixXd V = V0;
  MatrixXi F = F0;

  // declare input decimation
  Surface_mesh M_hat;
  MatrixXd V_coarse;
  MatrixXi F_coarse;

  // Loop over levels
  for(int i = 0;i<k;i++){
    // Use a provided input decimation if one was given for this level ...
    if(i < (int)V_decim.size() && V_decim[i].rows()>0)
    {
      V_coarse = V_decim[i];
      F_coarse = F_decim[i];
      // convert to CGAL format
      mesh_to_polyhedron(V_coarse,F_coarse,M_hat);
    }
    // ... otherwise compute it by decimating the previous layer with CGAL
    // decimation (regular or adaptive)
    else
    {
      const bool adaptive = !(i < (int)regular.size() && regular[i]);
      // specified number of faces for this cage
      const int Li = num_faces[i];
      // value to pass to the decimator
      float ratio = (1.*Li)/(1.*M.size_of_facets());
      // the previously computed cage will be decimated
      M_hat = M;
      // decimate
      decimate_CGAL(&M_hat,ratio,adaptive);
    }

    // Convert decimation to LibIGL/Eigen format
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
    // If input coarse mesh self-intersect, remove self-intersections with Meshfix
    if (IF.rows()>0)
    {
      #ifdef VERBOSE_DEBUG
        cout << i+1 << "-th input decimation self-intersects. Fixing with Meshfix " << endl;
      #endif
#if WITH_MESHFIX
      cout << "Polishing M" << i+1 << "..." << endl;
      meshfix_wrapper(MatrixXd(V_coarse),MatrixXi(F_coarse),V_coarse,F_coarse);
      cout << "Success!" << endl;
#else
      cout << "[WITH_MESHFIX not defined] Skipping polishing of M" << i+1 << "..." << endl;
#endif
      remesh_self_intersections(V_coarse,F_coarse,params,tempV,tempF,IF,J,IM);
      if (IF.rows()==0)
      {
        #ifdef VERBOSE_DEBUG
          cout << "Meshfix succesfully removed self-intersections" << endl;
        #endif
      }
      else{
          cout << "Wasn't able to remove all input self-intersections. Quitting..." << endl;
          return false;
      }
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
    cout << "Flowing M" << i << " inside M" << i+1 << "..." << endl;
    if (!flow_fine_inside_coarse(V,F,V_coarse,F_coarse,A_qv,H))
    {
      cout << "Flow failed to take fine mesh inside coarse mesh after 1000 iterations. Quitting" << endl;
      return false;
    }
    cout << "Success!" << endl;

    // Reinflate and output to cage to C
    MatrixXd C;
    cout << "Reinflating M" << i << ", pushing M" << i+1 << "..." << endl;
    reinflate(H,F,V_coarse,F_coarse,energy_expansion.c_str(),energy_final.c_str(),C);
    cout << "Success!" << endl;

    // sanity check: cage should never self-intersect at this stage
    remesh_self_intersections(C,F_coarse,params,tempV,tempF,IF,J,IM);
    if (IF.rows()>0){
      cout << i+1 << "-th output cage self-intersects. ERROR! Quitting...  " << endl;
      return false;
    }

    // sanity check: cage should never intersect input decimation
    intersect_other(C,F_coarse,V,F,true,IF);
    if (IF.rows()>0){
      cout << i+1 << "-th output cage intersect previous cage. ERROR! Quitting...  " << endl;
      return false;
    }

    // output cage is the input for the next level
    M.clear();
    mesh_to_polyhedron(C,F_coarse,M);
    V = C;
    F = F_coarse;

    // record the cage for this level
    Vout.push_back(C);
    Fout.push_back(F_coarse);

    // stream this cage to the caller as soon as it is ready
    if(per_level_callback)
    {
      per_level_callback(i, C, F_coarse);
    }
  }

  return true;
}
