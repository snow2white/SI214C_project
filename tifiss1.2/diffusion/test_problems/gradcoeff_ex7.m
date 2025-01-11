function [dcoeffdx,dcoeffdy] = specific_gradcoeff(x,y,nel)
%GRADCOEFF_EX7 computes the gradient of the coefficients a(x,y)
%   [dcoeffdx,dcoeffdy] = specific_gradcoeff(xy,evt);
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements
%   output
%          dcoeffdx   first component derivative
%          dcoeffdy   second component derivative
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  dcoeffdx = -x./2;
  dcoeffdy = -y./2;

end   % end function
