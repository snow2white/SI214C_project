function [intres_primal,intres_dual] = goafem_intres_p1_with_p1(xy,xl_s,yl_s,evt,ugal,zgal,subdivPar)
%GOAFEM_INTRES_P1_WITH_P1 computes elementwise interior residuals for both primal and dual solutions
%
% [intres_primal,intres_dual] = goafem_intres_p1_with_p1(xy,xl_s,yl_s,evt,ugal,zgal,subdivPar)
%
% input:
%               xy    vertex coordinate vector
%        xl_s,yl_s    local coordinats of sub-elements
%              evt    element mapping matrix
%             ugal    primal P1 solution vector
%             zgal    dual P1 solution vector
%        subdivPar    red/bisec3 uniform sub-division switch
%
% output: 
%    intres_primal    interior residuals (primal problem)
%      intres_dual    interior residuals (dual problem)
%
% The function computes the elementwise 'internal' residuals of local 
% residual problems. For the primal solution this is given by
%
% intres_f := \int_K ( f0(x) + div(\vec{f}) ) v(x) dx 
%           + \int_K div( a(x)\grad ugal(x) ) v(x) dx 
%
% For the dual solution, in the expression above there would be g0, \vec{g} and 
% zgal, respectively.
%
% Function(s) called: triangular_gausspoints
%                     subelt_transf
%                     tderiv
%                     goafem_tgauss_gradcoeff
%                     goafem_tgauss_L2rhs
%                     goafem_tgauss_divH1
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel = size(evt,1);

% Construct the integration rule (3/7/19/73 Gaussian points)
  nngpt = 7;
  [s,t,wt] = triangular_gausspoints(nngpt);
   
% Recover local coordinates and local solutions
  xl_v   = zeros(nel,3);
  yl_v   = zeros(nel,3);
  usol_v = zeros(nel,3); 
  zsol_v = zeros(nel,3);
  for ivtx = 1:3
      xl_v(:,ivtx) = xy(evt(:,ivtx),1);
      yl_v(:,ivtx) = xy(evt(:,ivtx),2);
      usol_v(:,ivtx) = ugal(evt(:,ivtx));
      zsol_v(:,ivtx) = zgal(evt(:,ivtx));
  end
  
% Allocate memory
  intres_primal = zeros(nel,3); 
  intres_dual   = zeros(nel,3);
  bdem     = zeros(nel,4,3,3);
  bde      = zeros(nel,3,3);
  fdem     = zeros(nel,4,3);
  gdem     = zeros(nel,4,3);
  fde      = zeros(nel,3);
  gde      = zeros(nel,3);
  xl_m     = zeros(nel,3);
  yl_m     = zeros(nel,3);
  
% Loop over sub-elements
  for subelt = 1:4
      % Recover local coordinates over sub-element
      for ivtx = 1:3
          xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
          yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
      end     
      % Loop over Gauss points
      for igpt = 1:nngpt         
          sigpt = s(igpt);
          tigpt = t(igpt);
          wght = wt(igpt);       
          [sigptloc,tigptloc] = subelt_transf(sigpt,tigpt,subelt,subdivPar);
          %
          % Evaluate derivatives
          [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigptloc,tigptloc,xl_v,yl_v);
          [jac_m,~,phi_m,~,~] = tderiv(sigpt,tigpt,xl_m,yl_m); 
          %
          % Gradient of the diffusion coefficients
          [diffx,diffy] = goafem_tgauss_gradcoeff(sigpt,tigpt,xl_m,yl_m);
          %
          % L2 parts of the sources F and G
          [frhs_m,grhs_m] = goafem_tgauss_L2rhs(sigpt,tigpt,xl_m,yl_m);
          %
          % Divergence of the H1 parts of the sources F and G
          [fdiv_m,gdiv_m] = goafem_tgauss_divH1(sigpt,tigpt,xl_m,yl_m);
          %
          % Loop over mid-edge hat functions
          for j = 1:3  
              % Compute div(a*grad)-contribution = grad(a)*grad(u_h): loop
              % over vertices old hat functions
              for i = 1:3  
                  bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffx(:) .* dphidx_v(:,i) .* phi_m(:,j) .* invjac_v(:) .* jac_m(:);
                  bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffy(:) .* dphidy_v(:,i) .* phi_m(:,j) .* invjac_v(:) .* jac_m(:);
              end
              %
              % Compute rhs-contribution from the source f
              fdem(:,subelt,j) = fdem(:,subelt,j) + wght * frhs_m(:) .* phi_m(:,j) .* jac_m(:); % L2-part
              fdem(:,subelt,j) = fdem(:,subelt,j) + wght * fdiv_m(:) .* phi_m(:,j) .* jac_m(:); % div-part
              %
              % Compute rhs-contribution from the source g
              gdem(:,subelt,j) = gdem(:,subelt,j) + wght * grhs_m(:) .* phi_m(:,j) .* jac_m(:); % L2-part
              gdem(:,subelt,j) = gdem(:,subelt,j) + wght * gdiv_m(:) .* phi_m(:,j) .* jac_m(:); % div-part       
          end
          % end mid-edge hat functions loop
      end
      % end Gauss points loop
  end
% end sub-elements loop

% -----------------------------------------------------------------------------
% Manual assembly of subelement contributions
% -----------------------------------------------------------------------------
  if subdivPar == 1
      %
      % Red sub-division: assembling
      % 
      % First edge
      bde(:,1,1) = bdem(:,2,3,1) + bdem(:,3,2,1) + bdem(:,4,1,1);
      bde(:,1,2) = bdem(:,2,3,2) + bdem(:,3,2,2) + bdem(:,4,1,2);
      bde(:,1,3) = bdem(:,2,3,3) + bdem(:,3,2,3) + bdem(:,4,1,3);
      fde(:,1)   = fdem(:,2,3)   + fdem(:,3,2)   + fdem(:,4,1);
      gde(:,1)   = gdem(:,2,3)   + gdem(:,3,2)   + gdem(:,4,1);
      % Second edge
      bde(:,2,1) = bdem(:,1,3,1) + bdem(:,3,1,1) + bdem(:,4,2,1);
      bde(:,2,2) = bdem(:,1,3,2) + bdem(:,3,1,2) + bdem(:,4,2,2);
      bde(:,2,3) = bdem(:,1,3,3) + bdem(:,3,1,3) + bdem(:,4,2,3);
      fde(:,2)   = fdem(:,1,3)   + fdem(:,3,1)   + fdem(:,4,2);
      gde(:,2)   = gdem(:,1,3)   + gdem(:,3,1)   + gdem(:,4,2);
      % Third edge
      bde(:,3,1) = bdem(:,1,2,1) + bdem(:,2,1,1) + bdem(:,4,3,1);
      bde(:,3,2) = bdem(:,1,2,2) + bdem(:,2,1,2) + bdem(:,4,3,2);
      bde(:,3,3) = bdem(:,1,2,3) + bdem(:,2,1,3) + bdem(:,4,3,3);
      fde(:,3)   = fdem(:,1,2)   + fdem(:,2,1)   + fdem(:,4,3);            
      gde(:,3)   = gdem(:,1,2)   + gdem(:,2,1)   + gdem(:,4,3);
  else
      %
      % Bisec3 sub-division: assembling
      % 
      % First edge
      bde(:,1,1) = bdem(:,3,2,1) + bdem(:,4,2,1);
      bde(:,1,2) = bdem(:,3,2,2) + bdem(:,4,2,2);
      bde(:,1,3) = bdem(:,3,2,3) + bdem(:,4,2,3);
      fde(:,1)   = fdem(:,3,2)   + fdem(:,4,2);
      gde(:,1)   = gdem(:,3,2)   + gdem(:,4,2);    
      % Second edge
      bde(:,2,1) = bdem(:,1,3,1) + bdem(:,2,1,1) + bdem(:,3,3,1) + bdem(:,4,1,1); 
      bde(:,2,2) = bdem(:,1,3,2) + bdem(:,2,1,2) + bdem(:,3,3,2) + bdem(:,4,1,2);
      bde(:,2,3) = bdem(:,1,3,3) + bdem(:,2,1,3) + bdem(:,3,3,3) + bdem(:,4,1,3);
      fde(:,2)   = fdem(:,1,3)   + fdem(:,2,1)   + fdem(:,3,3)   + fdem(:,4,1);
      gde(:,2)   = gdem(:,1,3)   + gdem(:,2,1)   + gdem(:,3,3)   + gdem(:,4,1);
      % Third edge
      bde(:,3,1) = bdem(:,1,2,1) + bdem(:,2,2,1);
      bde(:,3,2) = bdem(:,1,2,2) + bdem(:,2,2,2);
      bde(:,3,3) = bdem(:,1,2,3) + bdem(:,2,2,3);
      fde(:,3)   = fdem(:,1,2)   + fdem(:,2,2);
      gde(:,3)   = gdem(:,1,2)   + gdem(:,2,2);
  end
  
% -----------------------------------------------------------------------------  
% Assemble interior residuals from rhs-contributions
% -----------------------------------------------------------------------------
  for i = 1:3
      intres_primal(:,i) = intres_primal(:,i) + fde(:,i);
      intres_dual(:,i)   = intres_dual(:,i) + gde(:,i);
  end

% Assemble by adding the div(a*grad)-contribution times Galerkin solution
  for j = 1:3
      for k = 1:3
          intres_primal(:,j) = intres_primal(:,j) + bde(:,j,k).*usol_v(:,k);
          intres_dual(:,j)   = intres_dual(:,j) + bde(:,j,k).*zsol_v(:,k);
      end
  end
  
end % end function