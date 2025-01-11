function [dZ] = zeta_fd(p)
%ZETA_FD Computes the (approximated) signed distance function for Z-shaped domain
%
% The output is needed as an input for distmesh2d
%
% functions called: ell2_fd       (computing the sign. dist. for L domain)
%                   trapezium_fd  (computing the sign. dist. for Trapezoid)
%                   dunion        (distmesh original P.O. Persson function)
%
%   TIFISS function: AB; 07 June 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi
%
%
% NOTE
% -------------------------------------------------------------------------
% The (approximated) signed distance function for the Z-shaped domain is
% given by computing the minimum between the signed distance function for
% a classic L-shaped domain and a rectangle trapezoid in the bottom right
% part of the domain


  dL = ell2_fd(p);
  dT = trapezoid_fd(p);
  dZ = dunion(dL,dT);

end % end function