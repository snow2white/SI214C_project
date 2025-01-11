function [g] = goafem_specific_L2goal(x,y,nel)
%GOAFEM_ZERO_L2GOAL zero L2 part of the RHS of the dual problem
%
% [f] = goafem_specific_L2goal(x,y,nel)
%
% input:
%        x     x coordinate vector
%        y     y coordinate vector 
%      nel     number of elements  
%
% output: 
%        g     L2 part of the RHS (dual)
%
% See also GOAFEM_ZERO_L2RHS
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  g = zeros(nel,1);

end % end function