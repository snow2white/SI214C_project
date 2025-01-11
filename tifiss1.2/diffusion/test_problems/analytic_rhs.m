function f = specific_rhs(x,y,nel)
%analytic_rhs   analytic RHS forcing function
%   f = specific_rhs(x,y,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements  
%   IFISS function: DJS; 18 July 2006.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
f= 2*y.*(1-y) +2* x.*(1-x);
return