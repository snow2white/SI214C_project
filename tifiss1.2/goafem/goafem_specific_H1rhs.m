function [f1,f2] = goafem_specific_H1rhs(x,y,nel)
%GOAFEM_ZERO_H1RHS zero H1 part of RHS of the primal problem
%
% [f1,f2] = goafem_specific_H1rhs(x,y,nel)
%   
% input:
%          x    x coordinate vector
%          y    y coordinate vector
%        nel    number of elements
%
% output: 
%    [f1,f2]    source vector of H1 part of rhs (primal)
%         
% See also GOAFEM_ZERO_H1GOAL
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  f1 = zeros(nel,1); 
  f2 = zeros(nel,1); 

end % end function