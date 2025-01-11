function [edgeres_primal,edgeres_dual] = goafem_edgeres_p1_with_p1(xy,evt,eboundt,ugal,zgal,eex,tve,els)
%GOAFEM_EDGERES_P1_WITH_P1 computes edge residuals for both primal and dual solutions
%
% [edgeres_primal,edgeres_dual] = goafem_edgeres_p1_with_p1(xy,evt,eboundt,ugal,zgal,eex,tve,els)
%
% input:
%               xy    vertex coordinate vector
%              evt    element mapping matrix
%          eboundt    element boundary mapping matrix
%             ugal    primal P1 solution vector
%             zgal    dual P1 solution vector            
%              eex    element connectivity array
%              tve    edge location array
%              els    elementwise edge lengths
%
% output:
%   edgeres_primal   edge residuals (primal problem)
%     edgeres_dual   edge residuals (dual problem)
%
% The function computes the 'edge' residuals of local residual problems. 
% For the primal solution this is given by
% 
% edgeres := (1/2)\int_{\partial K}(Jump(\vec{f}) + a(s)Jump(ugal))v(s) ds
%
% For the dual solution, in the expression above there would be \vec{g} and zgal, 
% respectively.
%
% Function(s) called: gausspoints_oned
%                     tderiv          
%                     goafem_p1fluxjmps
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel = size(evt,1);
  
% Construct the integration rule (1/2/3/4/7/10 Gaussian points)
  ngpt = 7;
  [oneg,onew] = gausspoints_oned(ngpt);
  
% Recover local coordinates  
  xl_v   = zeros(nel,3);
  yl_v   = zeros(nel,3);
  usol_v = zeros(nel,3);
  zsol_v = zeros(nel,3);
  for ivtx = 1:3
      xl_v(:,ivtx)   = xy(evt(:,ivtx),1);
      yl_v(:,ivtx)   = xy(evt(:,ivtx),2);
      usol_v(:,ivtx) = ugal(evt(:,ivtx));
      zsol_v(:,ivtx) = zgal(evt(:,ivtx));
  end
  
% Allocate memory
  edgeres_primal = zeros(nel,3);
  edgeres_dual   = zeros(nel,3);

% Recover first the 3 external normals of each element
  [nx,ny] = get_normals(xl_v,yl_v,els);
  	  
% Edge residual due to finite element solution
  for igpt = 1:ngpt
      sigpt = oneg(igpt);
      sigpt_ref_edge = (1.0 + sigpt)/2.0; % [-1,1] -> [0,1]    reference-edge map
      sigpt_l = (1.0 + sigpt)/4.0;        % [-1,1] -> [0,1/2]   LEFT sub-edge map
      sigpt_r = (3.0 + sigpt)/4.0;        % [-1,1] -> [1/2,1]  RIGHT sub-edge map
      wigpt = onew(igpt);
      %   
      % First edge
      [~,~,phi_v_1,~,~] = tderiv(sigpt_ref_edge,1-sigpt_ref_edge,xl_v,yl_v);
      %
      % Second edge
      [~,~,phi_v_2,~,~] = tderiv(0,sigpt_ref_edge,xl_v,yl_v);
      %
      % Third edge
      [~,~,phi_v_3,~,~] = tderiv(sigpt_ref_edge,0,xl_v,yl_v);
      %        
      % Jump of primal Galerkin solution (a(s)Jump(ugal)) over the LEFT and RIGHT sub-element edges
      [njmp_l_primal] = goafem_p1fluxjmps(xl_v,yl_v,usol_v,nx,ny,eboundt,evt,eex,tve,sigpt_l);         
      [njmp_r_primal] = goafem_p1fluxjmps(xl_v,yl_v,usol_v,nx,ny,eboundt,evt,eex,tve,sigpt_r);
      %
      % Jump of dual Galerkin solution (a(s)Jump(zgal)) over the LEFT and RIGHT sub-element edges
      [njmp_l_dual]   = goafem_p1fluxjmps(xl_v,yl_v,zsol_v,nx,ny,eboundt,evt,eex,tve,sigpt_l);         
      [njmp_r_dual]   = goafem_p1fluxjmps(xl_v,yl_v,zsol_v,nx,ny,eboundt,evt,eex,tve,sigpt_r);
      %
      % Assembling contributions:
      %
      % Contribution from the first edge
      edgeres_primal(:,1) = edgeres_primal(:,1) + (1/2) * wigpt * njmp_l_primal(:,1) .* phi_v_1(:,2) .* (els(:,1)/4);
      edgeres_primal(:,1) = edgeres_primal(:,1) + (1/2) * wigpt * njmp_r_primal(:,1) .* phi_v_1(:,3) .* (els(:,1)/4);
      edgeres_dual(:,1)   = edgeres_dual(:,1)   + (1/2) * wigpt * njmp_l_dual(:,1)   .* phi_v_1(:,2) .* (els(:,1)/4);
      edgeres_dual(:,1)   = edgeres_dual(:,1)   + (1/2) * wigpt * njmp_r_dual(:,1)   .* phi_v_1(:,3) .* (els(:,1)/4);
      %
      % Contribution from the second edge
      edgeres_primal(:,2) = edgeres_primal(:,2) + (1/2) * wigpt * njmp_l_primal(:,2) .* phi_v_2(:,3) .* (els(:,2)/4);
      edgeres_primal(:,2) = edgeres_primal(:,2) + (1/2) * wigpt * njmp_r_primal(:,2) .* phi_v_2(:,1) .* (els(:,2)/4);
      edgeres_dual(:,2)   = edgeres_dual(:,2)   + (1/2) * wigpt * njmp_l_dual(:,2)   .* phi_v_2(:,3) .* (els(:,2)/4);
      edgeres_dual(:,2)   = edgeres_dual(:,2)   + (1/2) * wigpt * njmp_r_dual(:,2)   .* phi_v_2(:,1) .* (els(:,2)/4);
      %
      % Contribution from the third edge
      edgeres_primal(:,3) = edgeres_primal(:,3) + (1/2) * wigpt * njmp_l_primal(:,3) .* phi_v_3(:,2) .* (els(:,3)/4);
      edgeres_primal(:,3) = edgeres_primal(:,3) + (1/2) * wigpt * njmp_r_primal(:,3) .* phi_v_3(:,1) .* (els(:,3)/4);
      edgeres_dual(:,3)   = edgeres_dual(:,3)   + (1/2) * wigpt * njmp_l_dual(:,3)   .* phi_v_3(:,2) .* (els(:,3)/4);
      edgeres_dual(:,3)   = edgeres_dual(:,3)   + (1/2) * wigpt * njmp_r_dual(:,3)   .* phi_v_3(:,1) .* (els(:,3)/4);
  end
% end loop over Gaussian points
   
% Edge residual due to jumps of source vectors (primal and dual) 
  [edgeres_fjump] = goafem_jump_H1(xl_v,yl_v,nx,ny,evt,eboundt,els,eex,tve,1);
  [edgeres_gjump] = goafem_jump_H1(xl_v,yl_v,nx,ny,evt,eboundt,els,eex,tve,2);
   
% Final edge residuals (primal and dual) 
  edgeres_primal = edgeres_primal + edgeres_fjump;
  edgeres_dual   = edgeres_dual   + edgeres_gjump;
  
end  % end function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [nx,ny] = get_normals(xv,yv,els)
%Recover the 3 external unit normals per each element in the mesh: nx and ny
%contain the x-component and the y-component of the normals respectively

  sxe1 = (xv(:,3) - xv(:,2))./els(:,1);  sye1 = (yv(:,3) - yv(:,2))./els(:,1);  % unit tangential components 1st edge
  sxe2 = (xv(:,1) - xv(:,3))./els(:,2);  sye2 = (yv(:,1) - yv(:,3))./els(:,2);  % unit tangential components 2nd edge 
  sxe3 = (xv(:,2) - xv(:,1))./els(:,3);  sye3 = (yv(:,2) - yv(:,1))./els(:,3);  % unit tangential components 3rd edge
	  
  nxe1 = sye1;    nye1 = -sxe1;                                                 % unit normal components 1st edge
  nxe2 = sye2;    nye2 = -sxe2;                                                 % unit normal components 2nd edge
  nxe3 = sye3;    nye3 = -sxe3;                                                 % unit normal components 3rd edge
	  
% Store the x-component and the y-component of the normals
  nx = [nxe1, nxe2, nxe3];   	  
  ny = [nye1, nye2, nye3];
  
end % end child function