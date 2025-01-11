function [dT] = trapezoid_fd(p)
%TRAPEZOID_FD computes the (approximated) signed distance function for a rectangle trapezoid
%located in the bottom part of the square [-1,1]^2.
%
% This function is called by the function ZETA_FD that computes
% the signed distance function for Z-shaped domain
%
%   TIFISS function: AB; 07 June 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi
%
%
% NOTE: diagonal line's slope
% -------------------------------------------------------------------------
%  slope >> 1 gives right-bottom square   (total L-shaped domain  with the upper rectangle)
%  slope == 1 gives rectangle trapezoid   (total pure Zeta domain  "    "    "       "    )
%  slope < 1  gives trapezoid/rectangle   (total crack domain      "    "    "       "    )
%  slope << 1 gives bottom rectangle      (total square domain     "    "    "       "    )

  
  global slope 
  
% parameters
  xt = 0.5;    yt = -0.5;
  e = 1;       f = 1; 

  one = max( -((p(:,1) - xt)/3 + e/2) , p(:,1) - xt - e/2 );
  two = max( -((p(:,2) - yt)   + f/2) , p(:,2) - yt - f/2 );
  
  % diagonal line according to slope's value
  line = -slope*p(:,1) + p(:,2);
  dT = max( max(one,two) , line );
 
end  % end function