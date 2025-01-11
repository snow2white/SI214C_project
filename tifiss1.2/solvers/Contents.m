% SOLVERS
%
% Files
%   a_cdt               - matrix-vector product for scalar operator
%   amg_coarsen_plot    - function for visualisation of amg coarsening
%   amg_grids_setup     - performs algebraic multigrid setup phase
%   amg_smoother        - performs smoothing 
%   amg_smoother_params - generates structure for AMG smoother parameters
%   amg_smoother_setup  - generates smoother data for AMG
%   amg_v_cycle         - performs one AMG V-cycle
%   bsolve              - solves Galerkin system using backslash
%   helpme_it           - iterative solvers interactive help
%   ilu0                - incomplete factorization with no fill in
%   it_solve            - iterative solution of predefined steady problem
%   m_amgzt             - AMG preconditioner for scalar operator
%   m_amgzz             - AMG preconditioner with multiple v-cycles
%   m_diagt             - action of diagonal preconditioning operator
%   m_ilut              - incomplete LU preconditioner
%   m_masscheb          - mass matrix Chebyshev preconditioning operator
%   m_massdiag          - mass matrix diagonal preconditioning operator
%   m_nonet             - "no preconditioning" operator
%   milu0               - modified incomplete factorization with no fill-in
%   resplot             - plot residuals computed by iterative solvers
