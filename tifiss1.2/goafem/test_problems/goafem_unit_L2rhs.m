function [f] = goafem_specific_L2rhs(x,y,nel)
%GOAFEM_UNIT_L2RHS unit L2 part of the RHS of the primal problem
%
% [f] = goafem_specific_L2rhs(x,y,nel)
%
% input:
%        x     x coordinate vector
%        y     y coordinate vector 
%      nel     number of elements  
%
% output: 
%        f     L2 part of the RHS (primal)
%
% See also GOAFEM_UNIT_L2GOAL
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  f = ones(nel,1);

end % end function