function [dcoeffdx,dcoeffdy] = specific_gradcoeff(x,y,nel)
%ZEROS_GRADCOEFF returns the zero gradient for constant coefficients
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

  dcoeffdx = zeros(nel,1);  
  dcoeffdy = zeros(nel,1);
  
end   % end function
