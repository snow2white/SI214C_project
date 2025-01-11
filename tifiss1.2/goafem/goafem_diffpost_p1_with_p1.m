function [elerr_primal, err_est_primal, ...
          elerr_dual,   err_est_dual] = ...
          goafem_diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,ugal,zgal,evtY,xyY,subdivPar)
%GOAFEM_DIFFPOST_P1_WITH_P1 computes hierarchical YP error estimator for both primal and dual solutions
%
% [elerr_primal, err_est_primal, ...
%  elerr_dual,   err_est_dual] = ...
%  goafem_diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,ugal,zgal,evtY,xyY,subdivPar)
%
% input:
%               xy    vertex coordinate vector  
%              evt    element mapping matrix
%              eex    element connectivity array
%              tve    edge location array
%              els    elementwise edge lengths
%          eboundt    element boundary mapping matrix
%             ugal    primal P1 solution vector
%             zgal    dual P1 solution vector
%             evtY    element mapping matrix for midpoints
%              xyY    vertex coordinate vector for midpoints
%        subdivPar    red or bisec3 uniform sub-division flag
%
% output:
%     elerr_primal    vector of element indicators (primal)
%   err_est_primal    global error estimate        (primal)
%       elerr_dual    vector of element indicators (dual)
%     err_est_dual    global error estimate        (dual)
%
% The function solves the discrete problems for the eY estimators:
%
%   B0(eY,v) = F(v) - B(ugal,v),     for all v \in Y,          (1)
%   B0(eY,v) = G(v) - B(zgal,v),     for all v \in Y.          (2)
%
% using a standard element residual technique to construct the corresponding 
% local residual problems over each element K in the mesh; see [BS16, eq.(5.3)]
%
% Recall that the source F(v) (resp. G(v)) of (1) (resp. of (2)) follows the 
% representation of [MS09] (and [FPZ16]); see also GOAFEM_FEMP1_ADIFF.
%
% References:
%
% [BS16] Bespalov, Silvester, Efficient adaptive stochastic Galerkin methods for 
% parametric operator equations, SIAM J. Sci. Comput., 38(4)A2118-A2140, 2016;
%
% [MS09] Mommer, Stevenson, A goal-oriented finite element method with
% convergence rates, SIAM J. Numer. Anal., 47(2)861-866, 2009;
%
% [FPZ16] Feischl, Praetorius, van der Zee, An abstract analysis of optimal 
% goal-oriented adaptivity, SIAM J. Numer. Anal., 54(3)1423-1448, 2016.
%
% Function(s) called:  triangular_gausspoints
%                      tderiv
%                      goafem_gauss_coeff
%                      goafem_intres_p1_with_p1
%                      goafem_edgeres_p1_with_p1
%                      goafem_diffpost_p1_yp_bc
% 
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel = size(evt,1);
  
% Recover local coordinates
  xl_v = zeros(nel,3);
  yl_v = zeros(nel,3);
  for ivtx = 1:3
      xl_v(:,ivtx) = xy( evt(:,ivtx),1 );
      yl_v(:,ivtx) = xy( evt(:,ivtx),2 ); 
  end

% Construct the integration rule (3/7/19/73 Gaussian points)
  nngpt = 7;
  [s,t,wt] = triangular_gausspoints(nngpt);
  
% Allocate memory: number of linear hat functions per element=3
  adem = zeros(nel,4,3,3);
  ae   = zeros(nel,3,3);      
  xl_s = zeros(nel,4,3); 
  yl_s = zeros(nel,4,3);
  xl_m = zeros(nel,3);
  yl_m = zeros(nel,3);
 
% -------------------------------------------------------------------------  
% STEP 1: coordinates of midpoints: four sub-elements
% -------------------------------------------------------------------------
% First physical mid-edge point
  xedge1(:,1)= 0.5*(xl_v(:,2) + xl_v(:,3));    
  yedge1(:,1)= 0.5*(yl_v(:,2) + yl_v(:,3));
% Second physical mid-edge point
  xedge2(:,1)= 0.5*(xl_v(:,1) + xl_v(:,3));   
  yedge2(:,1)= 0.5*(yl_v(:,1) + yl_v(:,3));
% Third physical mid-edge point
  xedge3(:,1)= 0.5*(xl_v(:,1) + xl_v(:,2));   
  yedge3(:,1)= 0.5*(yl_v(:,1) + yl_v(:,2)); 
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
% STEP 2: contributions from subelements of LHS of problems (1)-(2)
% -----------------------------------------------------------------------------   
  for subelt = 1:4    
      % Recover local coordinates of sub-elements
      for ivtx = 1:3
          xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
          yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
      end
      % Loop over Gauss points
      for igpt = 1:nngpt     
          sigpt = s(igpt);
          tigpt = t(igpt);
          wght = wt(igpt);        
          % Evaluate derivatives
          [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigpt,tigpt,xl_m,yl_m);
          %
          % Evaluate diffusion coefficients
          [diffx,diffy] = goafem_tgauss_coeff(sigpt,tigpt,xl_m,yl_m);
          %
          % Loop over the three mid-edge linear functions
          for j = 1:3
              for i = 1:3                
                  adem(:,subelt,i,j) = adem(:,subelt,i,j) + wght * diffx(:) .* dphidx_v(:,i) .* dphidx_v(:,j) .* invjac_v(:);
                  adem(:,subelt,i,j) = adem(:,subelt,i,j) + wght * diffy(:) .* dphidy_v(:,i) .* dphidy_v(:,j) .* invjac_v(:);
              end
          end
          % end mid-edge linear functions loop    
      end
      % end of Gauss point loop  
  end
% end subdivided element loop 

% -----------------------------------------------------------------------------
% STEP 3: manual assembly of subelement contributions
% -----------------------------------------------------------------------------
  if subdivPar == 1
      %
      % Red sub-division
      % 
      % First edge
      ae(:,1,1) = adem(:,2,3,3) + adem(:,3,2,2) + adem(:,4,1,1);
      ae(:,1,2) = adem(:,3,2,1) + adem(:,4,1,2);
      ae(:,1,3) = adem(:,2,3,1) + adem(:,4,1,3);
      % Second edge
      ae(:,2,1) = adem(:,3,1,2) + adem(:,4,2,1);
      ae(:,2,2) = adem(:,1,3,3) + adem(:,3,1,1) + adem(:,4,2,2);
      ae(:,2,3) = adem(:,1,3,2) + adem(:,4,2,3);  
      % Third edge     
      ae(:,3,1) = adem(:,2,1,3) + adem(:,4,3,1);
      ae(:,3,2) = adem(:,1,2,3) + adem(:,4,3,2);
      ae(:,3,3) = adem(:,1,2,2) + adem(:,2,1,1) + adem(:,4,3,3);  
  else
      %
      % Bisec3 sub-division
      % 
      % First edge
      ae(:,1,1) = adem(:,3,2,2) + adem(:,4,2,2);
      ae(:,1,2) = adem(:,3,2,3) + adem(:,4,2,1);
      % ae(:,1,3) = empty
      % Second edge
      ae(:,2,1) = adem(:,3,3,2) + adem(:,4,1,2);
      ae(:,2,2) = adem(:,1,3,3) + adem(:,2,1,1) + adem(:,3,3,3) + adem(:,4,1,1);
      ae(:,2,3) = adem(:,1,3,2) + adem(:,2,1,2);  
      % Third edge     
      % ae(:,3,1) = empty 
      ae(:,3,2) = adem(:,1,2,3) + adem(:,2,2,1);
      ae(:,3,3) = adem(:,1,2,2) + adem(:,2,2,2);
  end
 
% ----------------------------------------------------------------------------- 
% STEP 4: RHS side of problems (1)-(2)
% -----------------------------------------------------------------------------   
% Local rhs of the local residual problems are given by the difference of
%
% res_int := \int_K(f0(x) + div(\vec{f}))v(x)dx + \int_K div(a(x)\grad ugal(x))v(x)dx 
%
% which is the 'internal residual' and
%
% res_edge := (1/2)\int_{\partial K}(Jump(\vec{f}) + a(s)Jump(ugal))v(s) ds
%
% which is the 'edge residual' for each element K in the mesh. 
% The same for the dual error problem.

% Element residuals (primal and dual)
  [intres_primal, intres_dual] = goafem_intres_p1_with_p1(xy,xl_s,yl_s,evt,ugal,zgal,subdivPar);  

% Edge residual: this is independent from the subdivision chosen (primal and dual)
  [edgeres_primal, edgeres_dual] = goafem_edgeres_p1_with_p1(xy,evt,eboundt,ugal,zgal,eex,tve,els);
  
  fprintf('Primal problem: intres = %7.4e;  edge_res = %7.4e\n',norm(intres_primal),norm(edgeres_primal));
  fprintf('Dual problem:   intres = %7.4e;  edge_res = %7.4e\n',norm(intres_dual),norm(edgeres_dual));

% Final rhs (primal and dual)
  rhs_primal = intres_primal - edgeres_primal;
  rhs_dual   = intres_dual   - edgeres_dual;  
 
% Impose Dirichlet boundary conditions 
  Ae = ae;
% ... on primal problem:
  [~,rhs_primal] = goafem_diffpost_p1_bc(Ae,rhs_primal,xy,evt,eboundt,evtY,xyY);
  
% ... on dual problem:
  [ae,rhs_dual]  = goafem_diffpost_p1_bc(Ae,rhs_dual,xy,evt,eboundt,evtY,xyY);
    
% NOTE that the lhs matrix ade is the same for both problems (hence it is
% returned only once)

% -----------------------------------------------------------------------------
% STEP 5: solving the system
% -----------------------------------------------------------------------------

% LDLT vectorised factorization
  [ae] = element_ldlt_factorization(ae);
  
% Solver forward-backward substitutions (primal and dual)
  solveprimal  = element_lusolve(ae,rhs_primal);
  solvedual    = element_lusolve(ae,rhs_dual);
  elerr_primal = solveprimal';
  elerr_dual   = solvedual';

% -----------------------------------------------------------------------------
% STEP 6: element indicators
% -----------------------------------------------------------------------------
% Elementwise indicators ||eYP||_B0^2 (primal and dual)
  elerr_sq_primal = zeros(nel,1);
  elerr_sq_dual   = zeros(nel,1);
  for ivtx = 1:3
      elerr_sq_primal(:) = elerr_sq_primal(:) + rhs_primal(:,ivtx) .* elerr_primal(ivtx,:)';
      elerr_sq_dual(:)   = elerr_sq_dual(:)   + rhs_dual(:,ivtx)   .* elerr_dual(ivtx,:)';
  end
% Elementwise indicators ||eYP||_B0 (primal and dual)
  elerr_primal = sqrt(elerr_sq_primal);
  elerr_dual   = sqrt(elerr_sq_dual);
  
% -----------------------------------------------------------------------------
% STEP 7: global error estimates
% -----------------------------------------------------------------------------  
% Global B0-norm of the estimators: ||eYP||_B0 (primal and dual)
  err_est_primal = norm(elerr_primal,2);
  err_est_dual   = norm(elerr_dual,2);
  
end  % end function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [ae] = element_ldlt_factorization(ae)
% Elementwise LDL^T factorization of lhs matrix
  nel   = size(ae,1);  % number of elements
  nnode = size(ae,2);  % Number of mid-edge linear functions (=3)
  nn    = nnode;
  dd    = zeros(nel,nn); 
  rr    = zeros(nel,nn);
  for kk=1:nn-1  
      for pp = 1:kk-1
          rr(1:nel,pp) = dd(1:nel,pp).*ae(1:nel,kk,pp);
      end
      %
      dd(1:nel,kk) = ae(1:nel,kk,kk);
      for pp = 1:kk-1
          dd(1:nel,kk)= dd(1:nel,kk) - ae(1:nel,kk,pp).*rr(1:nel,pp);
      end
      %
      for ii = kk+1:nn
          for pp = 1:kk-1
              ae(1:nel,ii,kk) = ae(1:nel,ii,kk) - ae(1:nel,ii,pp).*rr(1:nel,pp);
          end
          ae(1:nel,ii,kk) = ae(1:nel,ii,kk)./dd(1:nel,kk);
      end
  end
  %
  for pp = 1:nn-1
      rr(1:nel,pp) = dd(1:nel,pp).*ae(1:nel,nn,pp);
  end
  %
  dd(1:nel,nn) = ae(1:nel,nn,nn);
  %
  for pp = 1:nn-1
      dd(1:nel,nn) = dd(1:nel,nn) - ae(1:nel,nn,pp).*rr(1:nel,pp);
  end
  %
  % overwrite diagonal entries
  for kk=1:nn
      ae(1:nel,kk,kk) = dd(1:nel,kk);
  end

end % child function