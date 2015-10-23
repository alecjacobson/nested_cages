// BEGIN of libigl includes
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/boundary_facets.h>
#include <igl/circulation.h>
#include <igl/centroid.h>
#include <igl/slice_mask.h>
#include <igl/remove_unreferenced.h>
#include <igl/REDRUM.h>
#include <igl/matlab_format.h>
#include <igl/edges.h>
#include <igl/unique_edge_map.h>
#include <igl/unique.h>
#include <igl/writePLY.h>
#include <igl/pathinfo.h>
#include <igl/intersect.h>
#include <igl/list_to_matrix.h>
#include <igl/cgal/remesh_self_intersections.h>
#include <igl/cgal/polyhedron_to_mesh.h>
#include <igl/massmatrix.h>
#include <Eigen/Core>
#include <iostream>
#include <set>
#include <algorithm>
#include <cstdlib>
// END of libigl includes

// BEGIN of CGAL includes
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
// Extended polyhedron items which include an id() field
#include <CGAL/Polyhedron_items_with_id_3.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h> 
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point ;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Surface_mesh; 
typedef Surface_mesh::Halfedge_handle Halfedge_handle ;
typedef Surface_mesh::Vertex_handle   Vertex_handle ;
namespace SMS = CGAL::Surface_mesh_simplification ;
typedef SMS::Edge_profile<Surface_mesh> Profile ;   
// END of CGAL includes

// useful namespaces
using namespace std;
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

// The following is a Visitor that keeps track of the simplification process.
// In this example the progress is printed real-time and a few statistics are
// recorded (and printed in the end).
//
struct Stats
{
  Stats() 
    : collected(0)
    , processed(0)
    , collapsed(0)
    , non_collapsable(0)
    , cost_uncomputable(0) 
    , placement_uncomputable(0) 
  {} 
  
  std::size_t collected ;
  std::size_t processed ;
  std::size_t collapsed ;
  std::size_t non_collapsable ;
  std::size_t cost_uncomputable  ;
  std::size_t placement_uncomputable ; 
} ;

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh>
{
  My_visitor( Stats* s) : stats(s){} 

  // Called during the collecting phase for each edge collected.
  void OnCollected( Profile const&, boost::optional<double> const& )
  {
    ++ stats->collected ;
    // std::cerr << "\rEdges collected: " << stats->collected << std::flush ;
  }                
  
  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(Profile const&          
                 ,boost::optional<double> cost
                 ,std::size_t             initial
                 ,std::size_t             current
                 )
  {
    ++ stats->processed ;
    if ( !cost )
      ++ stats->cost_uncomputable ;
      
    // if ( current == initial )
    //   std::cerr << "\n" << std::flush ;
    // std::cerr << "\r" << current << std::flush ;
  }                
  
  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(Profile const&          
                   ,boost::optional<Point>  placement
                   )
  {
    if ( !placement )
      ++ stats->placement_uncomputable ;
  }                
  
  // Called for each edge which failed the so called link-condition,
  // that is, which cannot be collapsed because doing so would
  // turn the surface mesh into a non-manifold.
  void OnNonCollapsable( Profile const& )
  {
    ++ stats->non_collapsable;
  }                
  
  // Called AFTER each edge has been collapsed
  void OnCollapsed( Profile const&, Vertex_handle )
  {
    ++ stats->collapsed;
  }                
  
  Stats* stats ;
} ;

// taken from edge_collapse_enriched_polyhedron.cpp (CGAL)
void decimate_CGAL(Surface_mesh* surface_mesh, float ratio, bool adaptive){

  // Surface_mesh surface_mesh; 
  // std::ifstream is(filename) ; is >> surface_mesh ;

  // The items in this polyhedron have an "id()" field 
  // which the default index maps used in the algorithm
  // need to get the index of a vertex/edge.
  // However, the Polyhedron_3 class doesn't assign any value to
  // this id(), so we must do it here:
  int index = 0 ;
  
  for( Surface_mesh::Halfedge_iterator eb = (*surface_mesh).halfedges_begin()
     , ee = (*surface_mesh).halfedges_end()
     ; eb != ee
     ; ++ eb
     ) 
    eb->id() = index++;

  index = 0 ;
  for( Surface_mesh::Vertex_iterator vb = (*surface_mesh).vertices_begin()
     , ve = (*surface_mesh).vertices_end()
     ; vb != ve
     ; ++ vb
     )  
    vb->id() = index++;
    
  // Decimate to output 10% resolution mesh
    SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);
 
  Stats stats ;
  My_visitor vis(&stats) ;
    
  // The index maps are not explicitelty passed as in the previous
  // example because the surface mesh items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  if (adaptive){
    int r = SMS::edge_collapse
             (*surface_mesh
             ,stop
             ,CGAL::visitor      (vis)
             );
  } else {
  int r = SMS::edge_collapse
           (*surface_mesh
           ,stop
           ,CGAL::get_cost     (SMS::Edge_length_cost  <Surface_mesh>())
                 .get_placement(SMS::Midpoint_placement<Surface_mesh>())
                 .visitor      (vis)
           );
   }

  
  return;

}

int remove_all_chars_and_count(char* str, char c) {
    int removed = 0;
    char *pr = str, *pw = str;
    while (*pr) {
        *pw = *pr++;
        removed += (*pw == c);
        pw += (*pw != c);
    }
    *pw = '\0';
    return removed;
}

// function to check is a char * is an integer
bool legal_int(char *str) {
    while (*str)
        if (!isdigit(*str++))
            return false;
    return true;
}

int main(int argc, char * argv[])
{

  if (argc==1){
    cout<<"Usage: ./nested_cages [filename.(off)] L(1) L(2) ... L(k)"<<endl;
    cout<<"  L(1) L(2) ... L(k) is the number of faces for each cage"<<endl;
    return 0;
  }
  
  // number of layers
  int k = argc-3;
  cout << "number of layers = " << k << endl;

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


  // for now output CGAL decimations
  int L[k];
  bool adaptive = true;
  Surface_mesh M_hat;
  for(int i = 0;i<k;i++){
    std::ifstream is_file(argv[i+2]);
    // check if argv is a valid file
    if (is_file){
      is_file >> M_hat;
    } else{
      // first check if last charcater of argv[i+2] is r. If it is, drop the 'r' and adaptive = false
      adaptive = remove_all_chars_and_count(argv[i+2], 'r')==0;
      // else check if argv[i+2] is a valid integer (throw an error if it is not)
      if (!legal_int(argv[i+2])){
        cout << "you have to pass integer values or valid input deimatations"  << endl;
        cout << "the invalid argument you have passed is " << argv[i+2] << endl;
        return 0;
      }
      L[i] = atoi(argv[i+2]);
      float ratio = (1.*L[i])/(1.*M.size_of_facets()); 
      M_hat = M;
      decimate_CGAL(&M_hat,ratio,adaptive);

      // Check if decimations self-intersect. If they do, fix with Meshfix
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
      cout << "number of intersecting triangles " << IF.rows() << endl; 

      // If it does interesct, try to fix with Meshfix      

      // Flow M inside M_hat

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
