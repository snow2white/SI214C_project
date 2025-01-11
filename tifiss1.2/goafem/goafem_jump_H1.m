function [edgeresjmp] = goafem_jump_H1(xl_v,yl_v,nx,ny,evt,eboundt,els,eex,tve,type)
%GOAFEM_JUMP_H1 computes the jump components of the H1 part of the residual problem
%
% [edgeresjmp] = goafem_jump_H1(xl_v,yl_v,nx,ny,evt,eboundt,els,eex,tve,type)
%
% input:
%        xl_v,yl_v    local elementwise vertices coordinates
%            nx,ny    elementwise components of unit exterior normals
%              evt    element mapping matrix
%          eboundt    element boundary mapping matrix
%              els    elementwise edge lengths	
%              eex    element connectivity array
%              tve    edge location array
%             type    1 for primal, 2 for dual
%
% output:
%       edgeresjmp    jumps of source vector (\vec{f} or \vec{g}) 
%                     contributing to the edge residual
%
% Given the 'edge' residuals of local residual problems, i.e., for the primal 
% solution, 
%
% res_edge := (1/2)\int_{\partial K}(Jump(\vec{f}) + a(s)Jump(ugal))v(s) ds
%
% the function computes the contribution associated with the jump of the H1
% part of the source, i.e, 
%
% (1/2)\int_{\partial K} Jump(\vec{f}) v(s) ds.
%
% For the dual solution, in the expression above there would be \vec{g}. 
%
% The current implementation computes these contributions exactly without 
% numerical integration since we assume that both 
% \vec{f} = (f1,f2) (for primal problem) and 
% \vec{g} = (g1,g2) (for dual problem)  
% are costant (i.e, f1,f2,g1,g2 are constant). 
% In particular, this is the setup of Mommer-Stevenson(2009) example, where 
% f1 (resp. g1) is a characteristic function for which is not possible to 
% assign either 1 or 0 when gaussian points are located on edges along 
% the "singularity" line.
%
% Vectors \vec{f} and \vec{g} and computed below according to the problem.
%
% References: 
% [MS09] Mommer, Stevenson, A goal-oriented finite element method with
% convergence rates, SIAM J. Numer. Anal., 47(2)861-866, 2009;
%
% Function(s) called: goafem_specific_H1rhs
%                     goafem_specific_H1goal   
%
% See also GOAFEM_P1FLUXJMPS
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel = size(evt,1);

% -----------------------------------------------------------------------------
% Recover the H1 source vector  
% ----------------------------------------------------------------------------- 
% In Mommer-Stevenson(2009) example, source vectors (\vec{f} and \vec{g})
% contain characteristic functions \chi_T. They are computed by checking if an 
% internal point of each element belongs to the region T or not.
% Note that if vectors do not contains characteristic functions, the
% implementation still works fine.

% Internal middle points	
  xp = sum(xl_v,2)/3;
  yp = sum(yl_v,2)/3;
  
% Source vector \vec{f} (or \vec{g})  
  if type == 1
      [compvec1,compvec2] = goafem_specific_H1rhs(xp,yp,nel);   % Primal problem
  elseif type == 2
      [compvec1,compvec2] = goafem_specific_H1goal(xp,yp,nel);  % Dual problem
  end
  
% Source vector 
  ff = [compvec1,compvec2];
  
% -----------------------------------------------------------------------------      
% Flux and jump
% -----------------------------------------------------------------------------  
  flux = zeros(nel,4);
  jmp  = zeros(nel,3);
  
  flux(:,1) = ff(:,1) .* nx(:,1) + ff(:,2) .* ny(:,1);
  flux(:,2) = ff(:,1) .* nx(:,2) + ff(:,2) .* ny(:,2);
  flux(:,3) = ff(:,1) .* nx(:,3) + ff(:,2) .* ny(:,3);
  
% Replace zero indices in array tve by 4s
  tvx = tve;
  tvx(tvx==0) = 4;  
% Jump
  jmp(:,1) = flux(:,1) + flux( sub2ind([nel,4], eex(:,1), tvx(:,1)) );   % Jump on the first edge
  jmp(:,2) = flux(:,2) + flux( sub2ind([nel,4], eex(:,2), tvx(:,2)) );   % Jump of the second edge
  jmp(:,3) = flux(:,3) + flux( sub2ind([nel,4], eex(:,3), tvx(:,3)) );   % Jump of the third edge

% Remove Dirichlet boundary edge contributions: put 0 to the jump on boundary edges
  jmp( sub2ind([nel,3],eboundt(:,1),eboundt(:,2)) ) = 0.0;
  
% Since both components of sourcevec are constant, the jump will be also 
% constant and it goes outside the integral; then, it remains to compute the 
% integral \int_\Ej v_Ej(s) ds, j=1,2,3, with v_Ej(s) being the basis function 
% associated with the midpoint of the j-th edge. 
% NOTE that in P1 case, such integral is nothing but that the area of the 
% isosceles triangle with base=edge_j and height=1.   

  areae1 = 0.5 * els(:,1);  % \int_\E1 v_E1(s) ds = 0.5 * lenght(edge1) * 1
  areae2 = 0.5 * els(:,2);  % \int_\E1 v_E1(s) ds = 0.5 * lenght(edge2) * 1
  areae3 = 0.5 * els(:,3);  % \int_\E1 v_E1(s) ds = 0.5 * length(edge3) * 1
  
% Edge residual due to the source vector f
  edgeresjmp = zeros(nel,3);
  edgeresjmp(:,1) = (1/2) * jmp(:,1) .* areae1;
  edgeresjmp(:,2) = (1/2) * jmp(:,2) .* areae2;
  edgeresjmp(:,3) = (1/2) * jmp(:,3) .* areae3;
  
end  % end function