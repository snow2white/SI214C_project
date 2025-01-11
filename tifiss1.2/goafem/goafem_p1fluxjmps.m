function [jmp] = goafem_p1fluxjmps(xl_v,yl_v,sl_v,nx,ny,eboundt,evt,eex,tve,s)
%GOAFEM_P1FLUXJMPS vectorised edge jumps of P1 stochastic Galerkin solution
%
% [jmp] = goafem_p1fluxjmps(xl_v,yl_v,sl_v,nx,ny,eboundt,evt,eex,tve,s)
%	
% input:
%          xl_v       local elementwise x-coordinates
%	       yl_v       local elementwise y-coordinates 
%          sl_v       local elementwise P1 solution 
%            nx       x-component of exterior unit normals
%            ny       y-component of exterior unit normals
%       eboundt       element boundary mapping matrix
%           evt       element mapping matrix
%           eex       element connectivity array
%           tve       edge location array
%             s       gaussian point in [0,1]
% output:
%           jmp       component elementwise edge flux jumps
%
% Function(s) called: tderiv
%                     goafem_tgauss_coeff
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel = length(evt(:,1)); 

% -----------------------------------------------------------------------------
% Compute flux across edges
% -----------------------------------------------------------------------------

% Flux across the first 'reference' edge - diagonal
  [flux1] = compute_flux(s,1-s,xl_v,yl_v,sl_v,nx(:,1),ny(:,1),nel);
  
% Flux across the second 'reference' edge - left
  [flux2] = compute_flux(0,s,  xl_v,yl_v,sl_v,nx(:,2),ny(:,2),nel);  
		
% Flux across the third 'reference' edge - bottom
  [flux3] = compute_flux(s,0,  xl_v,yl_v,sl_v,nx(:,3),ny(:,3),nel);  
				  
% Allocate a zero columns for boundary jumps 
  flux = [flux1, flux2, flux3, zeros(nel,1)];
  
% -----------------------------------------------------------------------------  
% Compute jumps  
% ----------------------------------------------------------------------------- 
  jmp = zeros(nel,3); 

% Replace zero indices in array tve by 4s
  tvx = tve;
  tvx(tvx==0) = 4;
  jmp(:,1) = flux(:,1) + flux( sub2ind([nel,4],eex(:,1),tvx(:,1)) );  % Jump on the first edge
  jmp(:,2) = flux(:,2) + flux( sub2ind([nel,4],eex(:,2),tvx(:,2)) );  % Jump of the second edge
  jmp(:,3) = flux(:,3) + flux( sub2ind([nel,4],eex(:,3),tvx(:,3)) );  % Jump of the third edge

% Remove Dirichlet boundary edge contributions: put 0 to the jump on boundary edges	
  jmp( sub2ind([nel,3],eboundt(:,1),eboundt(:,2)) ) = 0.0;  

end  % end function
		
		
% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [flux] = compute_flux(xgp,ygp,xv,yv,sol,nxe,nye,nel)
%Computes the flux of the P1 finite element solution across the (1st, 2nd, or 3rd) edge;
%current gaussian point (xgp,ygp) is located on the (reference) edge.
  flux = zeros(nel,1);	  
  [~,invjac,~,dphidx,dphidy] = tderiv(xgp,ygp,xv,yv);
  [diffx,diffy] = goafem_tgauss_coeff(xgp,ygp,xv,yv);
  for j = 1:3
	  flux(:) = flux(:) + sol(:,j) .* ( diffx(:) .* dphidx(:,j) .* nxe + ...
	                                    diffy(:) .* dphidy(:,j) .* nye ) .* invjac(:);
  end 
  	  
end % end child function