function [A,F] = goafem_imposebc(A,F,xy,bound)
%GOAFEM_IMPOSEBC imposes Dirichlet boundary conditions
%
% [A,F] = goafem_imposebc(A,F,,xy,bound)
% 
% input:
%           A   stiffness matrix (without bcs imposed)
%           f   rhs vector (without bcs imposed)
%          xy   vertex coordinate vector  
%       bound   boundary vertex vector
%
% output:
%           A   stiffness matrix (with bcs imposed)
%           F   rhs vector (with bcs imposed)
%
% This function is based on original TIFISS function NONZEROBC 
% (DJS; 4 March 2005)
%
% Function(s) called: goafem_specific_bc
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nvtx     = length(F);
  nbd      = length(bound);
  null_col = sparse(nvtx,nbd);
  null_row = sparse(nbd,nvtx);
  
% Extract LHS matrix and RHS vector for interior nodes  
  Ax = A(1:nvtx,1:nvtx);
  fx = F(1:nvtx);
 
% Set boundary conditions
  xbd = xy(bound,1);
  ybd = xy(bound,2);
  bc  = goafem_specific_bc(xbd,ybd);

% Update RHS for interior nodes
  fx = fx - Ax(:,bound)*bc;
  fx(bound) = bc; 
  
% Update LHS for interior nodes
  dA          = zeros(nvtx,1);
  dA(bound)   = ones(nbd,1);
  Ax(:,bound) = null_col;
  Ax(bound,:) = null_row; 
  Ax = Ax + spdiags(dA,0,nvtx,nvtx);

% Update 
  A = Ax;
  F = fx;

end % end function