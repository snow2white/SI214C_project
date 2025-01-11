function [dL] = ell2_fd(p)
%ELL2_FD computes the (approximated) signed distance function for L-shaped domain
%
% This function is called by the function ZETA_FD that computes the signed
% distance function for Z-shaped domain.
% 
%   TIFISS function: AB; 07 June 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi
%
%
% NOTE:
% -------------------------------------------------------------------------
% The (approximated) distance function for the L-shaped domain is given by 
% the minimum of the distance functions for two rectangles, an upper
% horizontal rectangle (in [-1,1]x[0,1]) and a right vertical rectangle 
% (in [0,1]x[-1,1]), respectively.


  % Distance function for the upper horizontal rectangle      
  x0 = 0;        y0 = 0.5;                                    
  a = 2;         b = 1;                  
  ux1 = max( -((p(:,1)-x0)+a/2) , p(:,1)-x0-a/2);
  uy1 = max( -((p(:,2)-y0)+b/2) , p(:,2)-y0-b/2);
  uz1 = max(ux1,uy1);  
  
  % Distance function for the right vertical rectangle
  x01 = 0.5;     y01 = 0.0;  
  c = 1;         d = 2;  
  ux2 = max( -((p(:,1)-x01)+c/2) , p(:,1)-x01-c/2);
  uy2 = max( -((p(:,2)-y01)+d/2) , p(:,2)-y01-d/2);
  uz2 = max(ux2,uy2);

  % L-shaped domain distance
  dL = min(uz1,uz2);
 
end  % end function