function [divg] = goafem_specific_divH1goal(x,y,nel)
%GOAFEM_ZERO_DIVH1GOAL zero divergence of H1 part of RHS of the dual problem
%
% [divg] = goafem_specific_divH1rhs(x,y,nel)
%
% input:
%        x   x coordinate vector
%        y   y coordinate vector
%      nel   number of elements
%
% output: 
%     divg   divergence of source vector of the RHS (dual)
%
% The function returns the divergence of the source vector \vec{g} of the RHS 
% of the dual problem; see also GOAFEM_SPECIFIC_H1GOAL.
%
% See also GOAFEM_ZERO_DIVH1RHS
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  divg = zeros(nel,1);

end % end function