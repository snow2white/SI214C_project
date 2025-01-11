function [ade,bde,bdem,fde,gde] = goafem_diffpost_p1_with_p1_contrib(xy,evt,subdivPar)
%GOAFEM_DIFFPOST_P1_WITH_P1_CONTRIB elementwise contributions on both LHS and RHS of error problems
%
% [ade,bde,bdem,fde,gde] = goafem_diffpost_p1_with_p1_contrib(xy,evt,subdivPar)
%
% input:
%               xy     vertex coordinate vector
%              evt     element mapping matrix
%        subdivPar     red or bisec3 uniform sub-division flag
%
% output:
%              ade     elementwise LHS contribution 
%              bde     elementwise RHS contribution from Galerkin solutions
%             bdem     sub-elementwise RHS contribution from fe solution
%              fde     elementwise RHS contribution (primal)
%              gde     elementwise RHS contribution (dual)
%
% Function(s) called:  triangular_gausspoints
%                      tderiv
%                      goafem_tgauss_coeff
%                      goafem_tgauss_L2rhs
%                      goafem_tgauss_H1rhs
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel  = size(evt,1);   % Number of elements
  
% Recover local coordinates
  xl_v = zeros(nel,3); 
  yl_v = zeros(nel,3); 
  for ivtx = 1:3
      xl_v(:,ivtx) = xy( evt(:,ivtx), 1);
      yl_v(:,ivtx) = xy( evt(:,ivtx), 2);
  end
  
% Initialise local matrices  
  adem = zeros(nel,4,3,3);
  xl_s = zeros(nel,4,3); 
  yl_s = zeros(nel,4,3);
  xl_m = zeros(nel,3);
  yl_m = zeros(nel,3);
  
% ----------------------------------------------------------------------------- 
% STEP 1: coordinates of midpoints: four sub-elements
% -----------------------------------------------------------------------------
% First physical mid-edge point
  xedge1(:,1) = 0.5 * (xl_v(:,2) + xl_v(:,3));    
  yedge1(:,1) = 0.5 * (yl_v(:,2) + yl_v(:,3));
  
% Second physical mid-edge point
  xedge2(:,1) = 0.5 * (xl_v(:,3) + xl_v(:,1));  
  yedge2(:,1) = 0.5 * (yl_v(:,3) + yl_v(:,1));
  
% Third physical mid-edge point
  xedge3(:,1) = 0.5 * (xl_v(:,1) + xl_v(:,2));  
  yedge3(:,1) = 0.5 * (yl_v(:,1) + yl_v(:,2));

% Define the local sub-division 
  if subdivPar == 1
      %
      % Red sub-division
      % 
      % First physical sub-element 
      xl_s(:,1,1) = xl_v(:,1);      yl_s(:,1,1) = yl_v(:,1);
      xl_s(:,1,2) = xedge3(:);      yl_s(:,1,2) = yedge3(:);
      xl_s(:,1,3) = xedge2(:);      yl_s(:,1,3) = yedge2(:);
      % Second physical sub-element   
      xl_s(:,2,1) = xedge3(:);      yl_s(:,2,1) = yedge3(:);
      xl_s(:,2,2) = xl_v(:,2);      yl_s(:,2,2) = yl_v(:,2);
      xl_s(:,2,3) = xedge1(:);      yl_s(:,2,3) = yedge1(:);
      % Third physical sub-element 
      xl_s(:,3,1) = xedge2(:);      yl_s(:,3,1) = yedge2(:);
      xl_s(:,3,2) = xedge1(:);      yl_s(:,3,2) = yedge1(:);
      xl_s(:,3,3) = xl_v(:,3);      yl_s(:,3,3) = yl_v(:,3);
      % Fourth physical sub-element 
      xl_s(:,4,1) = xedge1(:);      yl_s(:,4,1) = yedge1(:);
      xl_s(:,4,2) = xedge2(:);      yl_s(:,4,2) = yedge2(:);
      xl_s(:,4,3) = xedge3(:);      yl_s(:,4,3) = yedge3(:);    
  else
      %
      % Bisec3 sub-division
      % 
      % First physical sub-element
      xl_s(:,1,1) = xl_v(:,1);      yl_s(:,1,1) = yl_v(:,1);
      xl_s(:,1,2) = xedge3(:);      yl_s(:,1,2) = yedge3(:);
      xl_s(:,1,3) = xedge2(:);      yl_s(:,1,3) = yedge2(:);
      % Second physical sub-element   
      xl_s(:,2,1) = xedge2(:);      yl_s(:,2,1) = yedge2(:);
      xl_s(:,2,2) = xedge3(:);      yl_s(:,2,2) = yedge3(:);
      xl_s(:,2,3) = xl_v(:,2);      yl_s(:,2,3) = yl_v(:,2);
      % Third physical sub-element 
      xl_s(:,3,1) = xl_v(:,2);      yl_s(:,3,1) = yl_v(:,2);
      xl_s(:,3,2) = xedge1(:);      yl_s(:,3,2) = yedge1(:);
      xl_s(:,3,3) = xedge2(:);      yl_s(:,3,3) = yedge2(:);
      % Fourth physical sub-element 
      xl_s(:,4,1) = xedge2(:);      yl_s(:,4,1) = yedge2(:);
      xl_s(:,4,2) = xedge1(:);      yl_s(:,4,2) = yedge1(:);
      xl_s(:,4,3) = xl_v(:,3);      yl_s(:,4,3) = yl_v(:,3);
  end    

% ----------------------------------------------------------------------------- 
% STEP 2: LHS of the linear system (the same for both primal and dual problems)
% ----------------------------------------------------------------------------- 

% Construct the integration rule
  nngpt = 7;
  [s,t,wt] = triangular_gausspoints(nngpt);

% Computing deterministic contribution from sub-elements 
  for subelt = 1:4 
      % Recover local coordinates of the current subelement
      for ivtx = 1:3
          xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
          yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
      end     
      % Loop over Gauss points
      for igpt = 1:nngpt         
          sigpt = s(igpt);
          tigpt = t(igpt);
          wght  = wt(igpt);
          % Evaluate derivatives and stochastic coefficients
          [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigpt,tigpt,xl_m,yl_m);  
          [diffx,diffy] = goafem_tgauss_coeff(sigpt,tigpt,xl_m,yl_m);
          % B0 bilinear form (i.e., with a_0(x) = coeff_v(:,1))
          for j = 1:3
              for i = 1:3                    
                  adem(:,subelt,i,j) = adem(:,subelt,i,j) + wght * diffx(:) .* dphidx_v(:,i) .* dphidx_v(:,j) .* invjac_v(:);
                  adem(:,subelt,i,j) = adem(:,subelt,i,j) + wght * diffy(:) .* dphidy_v(:,i) .* dphidy_v(:,j) .* invjac_v(:);              
              end
          end
      end
      % end of Gauss point loop
  end
% end of subdivided element loop

% Manual assembly of subelement contributions
  [ade] = assembling_lhs(adem,subdivPar);
  
% ----------------------------------------------------------------------------- 
% STEP 3: right-hand side of the linear system  
% ----------------------------------------------------------------------------- 
  
% Initialise local matrices
  bdem = zeros(nel,4,3,3); 
  fdem = zeros(nel,4,3);
  gdem = zeros(nel,4,3);
  xl_m = zeros(nel,3);
  yl_m = zeros(nel,3);
  
% Computing deterministic contributions from sub-elements  
  for subelt = 1:4
      for ivtx = 1:3
          xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
          yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
      end
      % Loop over Gauss points
      for igpt = 1:nngpt
          sigpt = s(igpt);   
          tigpt = t(igpt);
          wght  = wt(igpt);
          % Local transformation onto the sub-element of the reference element
          [sigptloc,tigptloc] = subelt_transf(sigpt,tigpt,subelt,subdivPar);
                  
          % Evaluate derivatives and coefficients
          [~,invjac_v,~,dphidx_v,dphidy_v]  = tderiv(sigptloc,tigptloc,xl_v,yl_v);
          [jac_m,~,phi_m,dphidx_m,dphidy_m] = tderiv(sigpt,tigpt,xl_m,yl_m);
          [diffx,diffy] = goafem_tgauss_coeff(sigpt,tigpt,xl_m,yl_m);

          % L2 rhs of both the primal and the dual problems
          [rhs_m,goal_m] = goafem_tgauss_L2rhs(sigpt,tigpt,xl_m,yl_m);

          % H1 rhs of both the primal and the dual problems
          [rhs1_m,rhs2_m,goal1_m,goal2_m] = goafem_tgauss_H1rhs(sigpt,tigpt,xl_m,yl_m);

          % Loop over sub-element vertices
          for j = 1:3  
              % Loop over X-basis functions
              % Contributions: \int_subelt a(x) \grad(Xbasis) \cdot \grad(Ybasis) dx 
              for i = 1:3                    
                  bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffx(:) .* dphidx_v(:,i) .* dphidx_m(:,j) .* invjac_v(:);                
                  bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffy(:) .* dphidy_v(:,i) .* dphidy_m(:,j) .* invjac_v(:);                
              end
              % Compute rhs-contributions:
              % Primal rhs
              fdem(:,subelt,j) = fdem(:,subelt,j) + wght * rhs_m(:)   .* phi_m(:,j) .* jac_m(:); % L2 part
              fdem(:,subelt,j) = fdem(:,subelt,j) - wght * rhs1_m(:)  .* dphidx_m(:,j);          % H1-part 1/2
              fdem(:,subelt,j) = fdem(:,subelt,j) - wght * rhs2_m(:)  .* dphidy_m(:,j);          % H1-part 2/2        
              % Dual rhs
              gdem(:,subelt,j) = gdem(:,subelt,j) + wght * goal_m(:)  .* phi_m(:,j) .* jac_m(:); % L2 part
              gdem(:,subelt,j) = gdem(:,subelt,j) - wght * goal1_m(:) .* dphidx_m(:,j);          % H1-part 1/2
              gdem(:,subelt,j) = gdem(:,subelt,j) - wght * goal2_m(:) .* dphidy_m(:,j);          % H1-part 2/2

          end
      end
      % end Gauss points loop
  end
% end sub-elements loop

% Manual assembly of sub-element contributions
  [bde,fde,gde] = assembling_rhs(bdem,fdem,gdem,subdivPar);
    
end  % end function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [ade] = assembling_lhs(adem,subdivPar)
% Elementwise assembling the contributions of the lhs for the 4 subelements

  ade  = zeros(size(adem,1),3,3);
  
  if subdivPar == 1
      % Red sub-division
      %
      % First edge
      ade(:,1,1) = adem(:,2,3,3) + adem(:,3,2,2) + adem(:,4,1,1);
      ade(:,1,2) = adem(:,3,2,1) + adem(:,4,1,2);
      ade(:,1,3) = adem(:,2,3,1) + adem(:,4,1,3);
      % Second edge
      ade(:,2,1) = adem(:,3,1,2) + adem(:,4,2,1);
      ade(:,2,2) = adem(:,1,3,3) + adem(:,3,1,1) + adem(:,4,2,2);
      ade(:,2,3) = adem(:,1,3,2) + adem(:,4,2,3);  
      % Third edge     
      ade(:,3,1) = adem(:,2,1,3) + adem(:,4,3,1);
      ade(:,3,2) = adem(:,1,2,3) + adem(:,4,3,2);
      ade(:,3,3) = adem(:,1,2,2) + adem(:,2,1,1) + adem(:,4,3,3);  
      
  else
      % Bisec3 sub-division: 
      %
      % First edge
      ade(:,1,1) = adem(:,3,2,2) + adem(:,4,2,2);
      ade(:,1,2) = adem(:,3,2,3) + adem(:,4,2,1);
      % ae(:,1,3) = empty
      % Second edge
      ade(:,2,1) = adem(:,3,3,2) + adem(:,4,1,2);
      ade(:,2,2) = adem(:,1,3,3) + adem(:,2,1,1) + adem(:,3,3,3) + adem(:,4,1,1);
      ade(:,2,3) = adem(:,1,3,2) + adem(:,2,1,2);
      % Third edge
      % ae(:,3,1) = empty 
      ade(:,3,2) = adem(:,1,2,3) + adem(:,2,2,1);
      ade(:,3,3) = adem(:,1,2,2) + adem(:,2,2,2);
  end

end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [bde,fde,gde] = assembling_rhs(bdem,fdem,gdem,subdivPar) 
% Elementwise assembling the contributions of the rhs for the 4 subelements

  bde = zeros(size(bdem,1),3,3);
  fde = zeros(size(bdem,1),3);   % integrals of the source for the primal rhs
  gde = zeros(size(bdem,1),3);   % integrals of the source for the dual rhs
  
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

end % end child function  