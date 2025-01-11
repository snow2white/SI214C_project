% GRIDS
%
% Files
%   adjust_unstruct_mesh    - renumbering nodes (i.e., edges) in evt for unstructured meshes  
%   circle_domain           - unit circle triangular grid generator
%   circle_obstacle_fd      - define circle obstacle geometry for DISTMESH
%   circle_obstacle_fh      - define circle obstacle density for DISTMESH
%   crack_domain_large      - large crack domain (crack on the left) structured grid generator
%   crack_domain_unit       - unit crack domain (crack on the right) structured grid generator
%   ell_domain              - L-shape structured grid generator
%   ell_domain_unstructured - L-shape unstructured grid generator
%   ell_fd                  - define the L-shape geometry for DISTMESH
%   ell_fh                  - set mesh density of L-shape geometry for DISTMESH
%   femlab_psfem            - convert FEMLAB mesh to TIFISS format
%   find_boundary_elements  - constructs eboundt and bound for Dirichlet BC
%   fint                    - subdivision function
%   fitint                  - computes contraction/expansion ratio
%   imeshplot               - triangular mesh verification
%   mesh_regularity_check   - ensures each element has at most 1 edge on boundary
%   natural_out_flowbc      - treats the right boundary nodes to be interior nodes
%   obstacle_domain         - obstacle domain triangular grid generator
%   p1_refinement           - uniformly refine the triangular mesh
%   p1grid                  - linear element grid generator
%   p1grid_detail_space     - linear detail space Y grid-generator
%   p2_grid_generator       - quadratic element grid generator
%   p2grid                  - quadratic element grid generator
%   p2p1grid                - P2-P1 element grid generator
%   plot_mesh_and_yspace    - plots mesh and its uniform (red/bisec3) refinement 
%   square_domain           - square domain Q2 grid generator
%   square_hole_domain      - square domain with central square hole structured grid generator 
%   square_obstacle_fd      - define square obstacle geometry for DISTMESH
%   square_obstacle_fh      - define square obstacle geometry for DISTMESH
%   subint                  - geometrically stretched subdivision generator
%   tedgegen                - edge information for flux jump computation
%   tedgegen_old            - legacy code
%   test_domain             - P1 grid generator for half L-shape
%   uniform_refinement      - uniformly refines the triangular mesh using either red or bisec3 refinement
