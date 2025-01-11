function f = specific_rhs(x,y,nel)
%RHS_EX6   RHS forcing function
%   f = specific_rhs(x,y,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements  
%   IFISS function: DJS; 28 February 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 

f = ( 2*ones(nel,1) - x.^2 - y.^2 )/8.0;

return
