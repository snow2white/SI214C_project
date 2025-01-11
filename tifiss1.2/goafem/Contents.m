% GOAFEM
%
% Files
%   goafem_convplot                    - plots the computed error estimates vs dofs
%   goafem_diffpost_p1_bc              - imposes Dirichlet BC at boundary elements for the residual problem
%   goafem_diffpost_p1_with_p1         - computes hierarchical YP error estimator for both primal and dual solutions
%   goafem_diffpost_p1_with_p1_2level  - two-level error estimation for both primal and dual solutions
%   goafem_diffpost_p1_with_p1_contrib - elementwise contributions on both LHS and RHS of error problems
%   goafem_diffpost_p1_with_p1_linsys  - computes spatial error estimation solving the fully assembled linear system
%   goafem_diffusion_main              - main driver for deterministic goal-oriented FEM
%   goafem_edgeres_p1_with_p1          - computes edge residuals for both primal and dual solutions
%   goafem_femp1_adiff                 - assembled P1 coefficient matrix generator for both primal and dual problems
%   goafem_imposebc                    - imposes Dirichlet boundary conditions
%   goafem_intres_p1_with_p1           - computes elementwise interior residuals for both primal and dual solutions
%   goafem_jump_H1                     - computes the jump components of the H1 part of the residual problem
%   goafem_marking_strategy            - Doerfler marking of elements/edges in goal-oriented framework 
%   goafem_mollifier                   - computes the normalization constant of [PO99] mollifier
%   goafem_p1fluxjmps                  - vectorised edge jumps of P1 stochastic Galerkin solution
%   goafem_plot_data                   - plots solution and spatial estimates
%   goafem_specific_bc                 - specification of Dirichlet boundary conditions
%   goafem_specific_coeff              - specification of diffusion operator 
%   goafem_specific_divH1goal          - zero divergence of H1 part of RHS of the dual problem
%   goafem_specific_divH1rhs           - zero divergence of H1 part of RHS of the primal problem
%   goafem_specific_gradcoeff          - zero gradient for constant diffusion operator
%   goafem_specific_H1goal             - Mommer-Stevenson(2009)-type H1 part of the RHS of the dual problem 
%   goafem_specific_H1rhs              - Mommer-Stevenson(2009)-type H1 part of the RHS of the primal problem 
%   goafem_specific_L2goal             - unit L2 part of the RHS of the dual problem
%   goafem_specific_L2rhs              - unit L2 part of the RHS of the primal problem
%   goafem_tgauss_coeff                - evaluates diffusion coefficient at triangle Gauss point
%   goafem_tgauss_divH1                - evaluates the divergence of the H1 part of the RHS of both primal and dual problems at Gaussian point
%   goafem_tgauss_gradcoeff            - evaluates derivatives of the coefficient at triangle Gauss point
%   goafem_tgauss_H1rhs                - evaluates the H1 part of the RHS of both the primal and dual problems at Gaussian point
%   goafem_tgauss_L2rhs                - evaluates the L2 part of the RHS of both the primal and dual problems at Gaussian point
%   triquad                            - Gaussian Quadrature for a triangular domain
