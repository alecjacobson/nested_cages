Matlab Repository for 'Nested Cages' project,
by Leonardo Sacht, Alec Jacobson and Etienne Vouga

1) Locate in 'usr/bin' the following binaries:
edge_collapse_enriched_polyhedron, SDFGen

2) This repository depends on other repositories
(maintained by Alec Jacobson): volume (winding number), gptoolbox

3) How to run in Matlab:

3.1) Read a mesh 'mesh.off' into the (V,F) format
>> [V0 F0] = readOFF('mesh.off');

3.2) Run 'multires_per_layer' specifying number of faces
per output level. In this example, we sepcify 3 levels halving the
number of faces from one level to the other:
>> [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'quadrature_order',2);

3.3) (cages_V{1}, cages_F{1}) is the coarsest level,
..., (cages_V{N}, cages_F{N}) is the finest level.
Now you can save these meshes to .OFF, visualize them
with meshplot, tetrahedralize them with Tetgen, etc.

3.4) For more information, issue
>> help multires_per_layer


## Installation

### Dependencies

 - gptoolbox
   - meshfix
 - cgal
 - eltopo
   - SuiteSparse
 - libigl (for c++ mex functions)


### Compile mex functions and command line programs

#### `edge_collapse_enriched_polyhedron`
This is a command line program (should be a mex, instead) created from the
CGAL examples. See top of `edge_collapse_enriched_polyhedron.cpp` for
compilation instructions.

(Might be out of date, see decimate_cgal)


#### `collide_eltopo_mex`
This is a mex function. 

Travel to `Code/eltopo/eltopo3d`. Build eltopo. You may need to edit
`Makefile.local_defs`. Then mex `collide_eltopo_mex.cpp`, see top of file for
instructions.

Then add path to make this function visible. Something like:

    addpath('../eltopo/eltopo3d/');

#### gptoolbox mex files
Set up matlab's mexopts with `-std=c++11 -stdlib=libc++` flags

Must compile:

    gptoolbox/mex/signed_distance.cpp
    gptoolbox/mex/intersect_other.cpp
    gptoolbox/mex/selfintersect.cpp
    gptoolbox/mex/in_element_aabb.cpp
    gptoolbox/mex/point_mesh_squared_distance.cpp
    gptoolbox/mex/winding_number.cpp
    gptoolbox/mex/decimate_cgal.cpp
