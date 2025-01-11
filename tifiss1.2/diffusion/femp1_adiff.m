function [a,q,f] = femp1_adiff(xy,evt)
%FEMP1_ADIFF set up linear variable diffusion matrices
%   [a,q,f] = femp1_adiff(xy,evt);
%   input
%          xy         vertex coordinate vector  
%          evt        element mapping matrix
%   output
%          a          diffusion matrix
%          q          mass matrix 
%          f          rhs vector
%    TIFISS function: DJS; 3 March 2017.
% Copyright (c) 2017 A. Bespalov, D.J. Silvester

% NOTE: difference with the original femp1_diff.m function
% -------------------------------------------------------------------------
% Here, the stiffness (and mass) matrices are computed in a loop over
% the set of gauss points because of the introduction of the
% non-constant coefficients (diffx,diffy)

  x = xy(:,1); 
  y = xy(:,2);
  nvtx = length(x);
  nel = size(evt,1);
  fprintf('setting up P1 diffusion matrices...  ');

% Initialise global matrices
  a = sparse(nvtx,nvtx);
  q = sparse(nvtx,nvtx);
  f = zeros(nvtx,1);

  xl_v = zeros(nel,3);
  yl_v = zeros(nel,3);
% Inner loop over elements    
  for ivtx = 1:3
      xl_v(:,ivtx) = x(evt(:,ivtx));
      yl_v(:,ivtx) = y(evt(:,ivtx)); 
  end
  
% Initialise local matrices
  ae = zeros(nel,3,3);
  qe = zeros(nel,3,3);
  fe = zeros(nel,3); 
   
% Construct the integration rule
  nngpt = 19;
  [s,t,wt] = triangular_gausspoints(nngpt);
  
% Loop over gaussian points  
  for igpt = 1:nngpt         
      sigpt = s(igpt);
      tigpt = t(igpt);
      wght = wt(igpt);       
    
      % Evaluate derivatives, rhs, and diffusion coefficients
      [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
      [rhs] = tgauss_source(sigpt,tigpt,xl_v,yl_v);
      [diffx,diffy] = tgauss_adiff(sigpt,tigpt,xl_v,yl_v);   
      
      for j = 1:3
          for i = 1:3
              % Stiffness matrix
              ae(:,i,j) = ae(:,i,j) + wght * diffx(:) .* dphidx(:,i) .* dphidx(:,j) .* invjac(:);
              ae(:,i,j) = ae(:,i,j) + wght * diffy(:) .* dphidy(:,i) .* dphidy(:,j) .* invjac(:);  
              % Mass matrix              
              qe(:,i,j) = qe(:,i,j) + wght * phi(:,i) .* phi(:,j) .* jac(:);
          end
          fe(:,j) = fe(:,j) + wght * rhs(:) .* phi(:,j) .* jac(:);
      end
      
  end

% Perform assembly of global matrix and source vector
  for krow = 1:3
      nrow = evt(:,krow);	 
      for kcol = 1:3
          ncol = evt(:,kcol);	  
          a = a + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
          q = q + sparse(nrow,ncol,qe(:,krow,kcol),nvtx,nvtx);
      end
%     f(nrow,1) = f(nrow,1) + fe(:,krow)
      for els = 1:nel
          f(nrow(els),1) = f(nrow(els),1) + fe(els,krow);
      end
  end
  
  fprintf('done\n')

end  % end function


