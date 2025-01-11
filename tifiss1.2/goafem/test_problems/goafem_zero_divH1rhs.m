function [divf] = goafem_specific_divH1rhs(x,y,nel)
%GOAFEM_ZERO_DIVH1RHS zero divergence of H1 part of RHS of the primal problem
%
% [divf] = goafem_specific_divH1rhs(x,y,nel)
%
% input:
%        x   x coordinate vector
%        y   y coordinate vector
%      nel   number of elements
%
% output: 
%     divf   divergence of source vector of the RHS (primal)
%
% The function returns the divergence of the source vector \vec{f} of the RHS 
% of the primal problem; see also GOAFEM_SPECIFIC_H1RHS.
%
% See also GOAFEM_ZERO_DIVH1GOAL
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  divf = zeros(nel,1);

end % end function