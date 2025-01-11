function [A,Q,F,G] = goafem_femp1_adiff(xy,evt)
%GOAFEM_FEMP1_ADIFF assembled P1 coefficient matrix generator for both primal and dual problems
%
% [A,Q,F,G] = goafem_femp1_adiff(xy,evt)
%
% input:
%         xy     vertex coordinate vector  
%        evt     element mapping matrix
%
% output:
%          A     cell structure with noarv matrices 
%          Q     mass matrix 
%          F     rhs vector for primal problem
%          G     rhs vector for dual problem
%
% This function extends FEMP1_ADIFF to the goal-oriented setting.
%
% The representation of the rhs (for both the primal and the dual problems)
% follows the one used in [MS09] (and [FPZ16]), which allows the estimation of 
% goal-functionals including (partial) derivatives of the (primal) solution 
% and also it introduces different non-geometric singularities in the primal 
% and/or dual solutions. 
%
% References:
%
% [MS09] Mommer, Stevenson, A goal-oriented finite element method with
% convergence rates, SIAM J. Numer. Anal., 47(2)861-866, 2009;
%
% [FPZ16] Feischl, Praetorius, van der Zee, An abstract analysis of optimal 
% goal-oriented adaptivity, SIAM J. Numer. Anal., 54(3)1423-1448, 2016;
% 
% Function(s) called: tderiv
%                     goafem_tgauss_coeff
%                     goafem_tgauss_L2rhs
%                     goafem_tgauss_H1rhs
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nvtx = size(xy,1);
  nel  = size(evt,1);
  x    = xy(:,1); 
  y    = xy(:,2);

% Initialise global matrices
  A = sparse(nvtx,nvtx);
  Q = sparse(nvtx,nvtx);
  F = zeros(nvtx,1);
  G = zeros(nvtx,1);
  
% Initialise local matrices
  ae = zeros(nel,3,3);
  qe = zeros(nel,3,3);
  fe = zeros(nel,3);
  ge = zeros(nel,3);

% Recover local coordinates  
  xl_v = zeros(nel,3);
  yl_v = zeros(nel,3);
  for ivtx = 1:3 
      xl_v(:,ivtx) = x(evt(:,ivtx));
      yl_v(:,ivtx) = y(evt(:,ivtx)); 
  end
  
% Construct the integration rule (3/7/19/73 Gaussian points)
  nngpt = 19;
  [s,t,wt] = triangular_gausspoints(nngpt);
  
% Loop over gaussian points  
  for igpt = 1:nngpt
      sigpt = s(igpt);
      tigpt = t(igpt);
      wght  = wt(igpt);
      % Evaluate derivatives
      [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
      %
      % Evaluate coefficients
      [diffx,diffy] = goafem_tgauss_coeff(sigpt,tigpt,xl_v,yl_v);
      %
      % Evaluate L2 part of the rhs
      [rhs,goal] = goafem_tgauss_L2rhs(sigpt,tigpt,xl_v,yl_v);
      % 
      % Evaluate H1 part of the rhs
      [rhs1,rhs2,goal1,goal2] = goafem_tgauss_H1rhs(sigpt,tigpt,xl_v,yl_v);
      %
      for j = 1:3
          for i = 1:3          
              % Stiffness matrix
              ae(:,i,j) = ae(:,i,j) + wght * diffx(:) .* dphidx(:,i) .* dphidx(:,j) .* invjac(:);
              ae(:,i,j) = ae(:,i,j) + wght * diffy(:) .* dphidy(:,i) .* dphidy(:,j) .* invjac(:);  
              % Mass matrix              
              qe(:,i,j) = qe(:,i,j) + wght * phi(:,i) .* phi(:,j) .* jac(:);        
          end
          % Primal rhs vector
          fe(:,j) = fe(:,j) + wght * rhs(:)   .* phi(:,j) .* jac(:); % L2-part
          fe(:,j) = fe(:,j) - wght * rhs1(:)  .* dphidx(:,j);        % H1-part 1/2
          fe(:,j) = fe(:,j) - wght * rhs2(:)  .* dphidy(:,j);        % H1-part 2/2
          % Dual rhs vector
          ge(:,j) = ge(:,j) + wght * goal(:)  .* phi(:,j) .* jac(:); % L2-part
          ge(:,j) = ge(:,j) - wght * goal1(:) .* dphidx(:,j);        % H1-part 1/2
          ge(:,j) = ge(:,j) - wght * goal2(:) .* dphidy(:,j);        % H1-part 2/2       
      end
  end

% Perform assembly of global matrix and source vectors
  for krow = 1:3
      nrow = evt(:,krow);	 
      for kcol = 1:3
          ncol = evt(:,kcol);	  
          A = A + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
          Q = Q + sparse(nrow,ncol,qe(:,krow,kcol),nvtx,nvtx);
      end
      for els = 1:nel
          F(nrow(els),1) = F(nrow(els),1) + fe(els,krow);
          G(nrow(els),1) = G(nrow(els),1) + ge(els,krow);
      end
  end
    
end % end function