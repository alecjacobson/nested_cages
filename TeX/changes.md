# list of changes made since submission

- cited and compares to [Sander et al.\ 2000] see \reffig{dane-vs-ben-chen}
- typos in text and equations noted by previous reviews
- handles.obj meshes have been fixed and checks have been added to
  algorithm to ensure that decimations meet assumptions on input
- clarified that multi-resolution Poisson equations are volumetric not
  surface-based and explicitly state that interior of cages are meshed with
  tetrahedra
- changed all \shortcite (e.g. ``Sander et al.\ ... [2000]") citations to \cite
  (e.g. ``[Sander et al.\ 2000] ...'')



# Changes that should be made with notes already in tex

 - Claim that watertightness is only important property need from decimation,
   but clear that homeomorphism is also needed 
 - Must check one vertex per connected component is inside
 - \cite{Baerentzen:2005:SDC} normals are not guaranteed to point inside for
   evaluation points _on_ the surface
 - In cases where coarse mesh is inflated, what keeps it tight? Energies in
   next section are defined w.r.t. original mesh
 - velocities/gradients in re-inflation should have minus signs?
 - Re-tet-meshing is needed for inhomogeneous equations where density is
   changing 
 - [Martin et al. 2008] used harmonic coordinates on general polyhedra for FEM.
 - overall timing performance should be better explained
 - cite \cite{Volino:2006:RSC} "Resolving Surface ...": could be used to replace
   re-inflation?
 - approximate gradient descent seems brute force: discuss Newton's method
   attempts: Only really makes sense for ARAP, doesn't make much of a
   difference and adds parameters (for stability).
 - I would have appreciated seeing (more) coarser cages
 - octopus in timings table differs from octopus used in multi-res example
 - no high genus inputs 
 - Cite papers that have used cell-replicating:
   \cite{Teran:2005:CSS,Nesme:2009:PTE}
 - previous multi-grid methods handle complicated boundaries
   \cite{Dick:2011,Ferstl:2014} Use adaptive octree multigrid solvers for
   fractures and fluids, respectively.  Though supporting complicated and even
   dynamic boundaries, they are not immune to troubles of regular grids:
   boundaries must be uniformly represented at the finest grid level to avoid
   aliasing.
 - Cite and discuss that surface offsetting and normal flowing is not a novel
   idea.
 - collision detection schemes do not need intermediary triangle so isosurface
   of distance fields are OK
 - Not fair to attack [Ben-Chen et al.] for non-convergence since ours might
   not
 - Frankenstein example complains that [Jacobson et al.] doesn't coarsen domain
   or create nice triangles, but problem statement in intro didn't include
   these goals
 - No discussion about output mesh quality: allude to ARAP energies in intro
 - No discussion of stopping criteria for re-inflation. Local minima? Closeness
   to global (feasible?) minimum? We could try running the method with fine
   mesh == coarse mesh.
 - no discussion of cons (pros are clear) of not incorporating decimation:
   losing extra combinatorial DOFs: mention experiments with subdivision,
   complicates formulation and probably enough to be its own contribution as
   future work
 - No reason for confidence that quadrature approximation of energy is "somehow
   as nice and useful as exact energy"
 - Flow is not guaranteed to even terminate. How is it stopped?
 - Doesn't discuss shrink wrapping alternative
 - Not clear why/that self-intersections of fine mesh during flow can be ignored
 - Give examples of "stupid choices" of energies that don't regularize flow
 - Octopi dropping doesn't say that they're rigid bodies
 - No discussion that nesting is not maintained during deformation of fine mesh
 - Simple example of naive decimation failing to converge for multi-resolution
 - Explain why 10 mins for backsubstitution? Are factors dense?
 - reference and compare to algebraic multigrid \cite{ruge1987algebraic}.
 - reference surface-based multi-resolution methods
   \cite{Aksoylu2005msu,Chuan:2009:ELO}
 - why not always use the point-to-plane distance as ICP community does?
 - Why are you comparing to conformalized mean curvature flow rather than
   standard MCF? cMCF is at least guaranteed to shrink to a point, MCF will
   converge to a zero-area, but not zero-length structure.
 - multiresolution results seem to suggest that without this method
   multiresolution is hopeless for the octopus and armadillo: it \emph{may} be
   possible to tweak extrapolation parameters to handle these cases with naive
   decimation but an automatic method for determining correct extrapolation for
   convergence is not obvious. We tried: linear extrapolation from nearest tet,
   constant interpolation of nearest vertex, linear interpolation of closest
   point on nearest face. Bottom line, extrapolation is not $C^0$, and some
   amount of cleverness is needed to avoid discontinuities.
 - what are the transfer operators used in multiresolution results? Appendix
   for multiresolution implementation? We use linear interpolation as
   prolongation operator and its transpose as the restriction operator (cite).
 - reference recent collision detection papers \cite{BEB2012,Wang:2014:DCC}: I
   guess we could at least say that we lean on implementations of collision
   detection in ElTopo/SPACM but we readily benefit from new advances in
   faster  robust detection.
 - Compare to direct solver on machine with more memory? No, there will always
   be a problem size big enough that direct solver will choke: This is just
   asymptotic analysis of memory.
 - How can point clouds/polygon soups be handled?

# Changes that should be made without notes in tex

 - Tone down homeomorphism requirement: topological noise
 - how does the method _rely_ on the input simplification?
 - when and in which configuration will the algorithm break?
 - "the initial flow is not guaranteed to work": show impossible csaszar-torus
   case, discuss types of meshes where flow fails: "surface too close to medial
   axis?"
 - theoretical snags are not adequately discussed
 - need stronger argument for strict nesting
 - frustrating that discussion of other options and previous work came after
   proposal of ideas
 - discuss drawbacks in situa rather than at the end
 - hahaha, "paper should be 8-pages instead of 12"
 - discussion of problem difficulty before jumping into concrete solution
