function [g1,g2] = goafem_specific_H1goal(x,y,nel)
%GOAFEM_ZERO_H1GOAL zero H1 part of RHS of the dual problem
%
% [g1,g2] = goafem_specific_H1rhs(x,y,nel)
%   
% input:
%          x    x coordinate vector
%          y    y coordinate vector
%        nel    number of elements
%
% output: 
%    [g1,g2]    source vector of H1 part of rhs (dual)
%         
% See also GOAFEM_ZERO_H1RHS
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  g1 = zeros(nel,1); 
  g2 = zeros(nel,1); 

end % end function