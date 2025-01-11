function [coeffx,coeffy] = goafem_specific_coeff(x,y,nel)
%GOAFEM_UNIT_COEFF unit constant diffusion operator 
%
% [coeffx,coeffy] = goafem_specific_coeff(x,y,nel)
%   
% input:
%                x   x coordinate vector
%                y   y coordinate vector 
%              nel   number of elements
%
% output:
%  [coeffx,coeffy]   diffusion coefficient
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  coeffx = ones(nel,1);
  coeffy = ones(nel,1); 

end % end function