#include "cgal.h"

// computes (overlapping) decimations using CGAL's routines
void decimate_CGAL(
  Surface_mesh* surface_mesh, 
  const float ratio, 
  const bool adaptive)

{

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
#if CGAL_VERSION_NR >= 1040700000
      ,CGAL::parameters::visitor      (vis)
#else
      ,CGAL::visitor      (vis)
#endif
      );
  } else {
  int r = SMS::edge_collapse
    (*surface_mesh
    ,stop
#if CGAL_VERSION_NR >= 1040700000
    ,CGAL::parameters::get_cost     (SMS::Edge_length_cost  <Surface_mesh>())
#else
    ,CGAL::get_cost     (SMS::Edge_length_cost  <Surface_mesh>())
#endif
    .get_placement(SMS::Midpoint_placement<Surface_mesh>())
    .visitor      (vis)
    );
   }

  
  return;

}
