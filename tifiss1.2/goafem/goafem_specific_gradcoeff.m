function [dcoeffdx,dcoeffdy] = goafem_specific_gradcoeff(x,y,nel)
%GOAFEM_ZERO_GRADCOEFF zero gradient for constant diffusion operator
%
% [dcoeffdx,dcoeffdy] = goafem_specific_gradcoeff(x,y,nel)
%   
% input:
%                   x   x coordinate vector
%                   y   y coordinate vector 
%                 nel   number of elements
%
% output:
% [dcoeffdx,dcoeffdy]   components derivative
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  dcoeffdx = zeros(nel,1);  
  dcoeffdy = zeros(nel,1);
  
end % end function