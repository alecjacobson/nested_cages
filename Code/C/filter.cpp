#include "filter.h"

// common
#include <array2.h>
#include <ccd_wrapper.h>
#include <collisionqueries.h>
#include <marching_tiles_hires.h>
#include <util.h>
#include <vec.h>
#include <wallclocktime.h>

// el topo
#include <collisionpipeline.h>
#include <eltopo.h>
#include <iomesh.h>
#include <meshrenderer.h>
#include <runstats.h>
#include <surftrack.h>
#include <trianglequality.h>

// talpa
#include <framestepper.h>
#include <meancurvature.h>
#include <meshdriver.h>
#include <scriptinit.h>
#include <simulation.h>

void reinflate(
  const Eigen::MatrixXd & Vf, 
  const Eigen::MatrixXi & T, 
  const Eigen::MatrixXd & Uf, 
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXi & F_hat, 
  Eigen::MatrixXd & Uc)
{
  using namespace Eigen;
  using namespace std;
  // Define concataned meshes (treat fine and coarse meshes as a big single mesh)
  MatrixXd V0,V1;
  V0.block(0,0,Vf.rows(),3) = Vf;
  V0.block(Vf.rows(),0,C.rows(),3) = C;
  V1.block(0,0,Vf.rows(),3) = Vf+Uf;
  V1.block(Vf.rows(),0,C.rows(),3) = C+Uc;
  MatrixXi F_all;
  F_all.block(0,0,T.rows(),3) = T;
  F_all.block(T.rows(),0,F_hat.rows(),3) = F_hat;

  // set vertex masses (=infty for fine mesh vertices 
  // and =1.0 for the last ones)
  double masses[V0.rows()];
  for (int i=0; i<Vf.rows(); i++){
      masses[i] = std::numeric_limits<double>::infinity();
  }
  for (int i=Vf.rows(); i<V0.rows(); i++){
      masses[i] = 1.0;
  }
  // encapsulate all data into an ElTopoMesh
  ElTopoMesh eltopo_time0;
  eltopo_time0.num_vertices = V0.rows();
  eltopo_time0.vertex_locations = V0.data();
  eltopo_time0.num_triangles = F_all.rows();
  eltopo_time0.triangles = F_all.data();
  eltopo_time0.vertex_masses = masses;

  double *V1a = V1.data();

  // Set general parameters
  ElTopoGeneralOptions sim_general_options;
  // do not print stuff to the console
  sim_general_options.m_verbose = 0;
  // do avoid self-intersections
  sim_general_options.m_collision_safety = 1;
  // separation between colliding meshes 
  sim_general_options.m_proximity_epsilon = 1e-3; 

  // Set Simulation parameters
  ElTopoIntegrationOptions sim_integration_options;
  sim_integration_options.m_friction_coefficient = 0.0;
  sim_integration_options.m_dt = 1.0;

  double* V_final;
  double out_dt = 0.0;
  el_topo_integrate(&eltopo_time0, V1a, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
  if (out_dt<1.0){
  	cout << "Failed to reach final positions out_dt = " << out_dt << endl;
  	return;
  }

  
}