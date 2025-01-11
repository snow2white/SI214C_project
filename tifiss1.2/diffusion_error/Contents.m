% DIFFUSION_ERROR
%
% Files
%   diffpost_p1_bc             - postprocesses Poisson error estimator at boundary elements
%   diffpost_p1_bc_x           - legacy code
%   diffpost_p1_with_p1        - a posteriori estimation for P1 using mid-edge P1 functions
%   diffpost_p1_with_p1_2level - computes 2-level error estimator for stochastic Galerkin P1 solution
%   diffpost_p1_with_p1_linsys - computes hierarchical eY estimator solving the (fully) assembled error problem
%   diffpost_p1_with_p2        - a posteriori estimation for P1 using P2 bubble functions
%   diffpost_p1_x              - legacy code
%   diffpost_p2_with_p4        - local Poisson error estimator for P2 with P4 correction
%   edgeres_p1_with_p1         - edge residuals for P1 solution using mid-edge P1 functions
%   edgeres_p1_with_p2         - edge residuals for P1 solution using P2 bubble functions
%   element_lusolve            - vectorized local backward-forward solves      
%   intres_p1_with_p1          - interior residuals for P1 solution using mid-edge P1 functions
%   intres_p1_with_p2          - interior residuals for P1 solution using P2 bubble functions
%   localbc_p                  - imposes Dirichlet BC for Poisson error estimator
%   localbc_p_new              - development code
%   p1fluxjmps                 - corrected flux jumps for triangular P1 grid
%   p1fluxjmps_x               - legacy code
%   p1res_diff                 - interior residuals for triangular P1 grid
%   p2fluxjmps                 - vectorised flux jumps for triangular P2 grid
%   reorder_s                  - locate 1D Gauss points on the reference element edges
%   reorder_s_x                - legacy code
%   subelt_transf              - locate 2D Gaussian points on the right reference subelement
%   tdiffpost_res              - computes P1 element residual error estimator
%   vtderiv                    - derivatives of linear shape functions for vector s, t
%   vtgauss_adiff              - elementwise permeability field for vector s, t
%   vtqarderiv                 - derivatives of quartic shape functions for vector s, t
%   vtqarshape                 - quartic shape functions for vector s, t
%   vtqderiv                   - derivatives of super-quadratic shape functions for vector s,t
%   vtqderiv_2                 - second derivatives of quadratic shape functions for vector s, t
%   vtqshape                   - super-quadratic shape functions for vector s, t
%   vtqshape_2                 - second derivatives on reference element for vector s, t
%   vtshape                    - linear shape functions for vector s, t
