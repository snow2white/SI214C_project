function [g] = goafem_specific_L2goal(x,y,nel)
%GOAFEM_UNIT_L2GOAL unit L2 part of the RHS of the dual problem
%
% [g] = goafem_specific_L2goal(x,y,nel)
%
% input:
%        x     x coordinate vector
%        y     y coordinate vector 
%      nel     number of elements  
%
% output: 
%        g     L2 part of the RHS (dual)
%
% See also GOAFEM_UNIT_L2GOAL
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  g = ones(nel,1);

end % end function