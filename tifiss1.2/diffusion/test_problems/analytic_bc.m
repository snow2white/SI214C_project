function bc = specific_bc(xbd,ybd)
%analytic_bc   Reference problem  boundary condition 
%   bc = specific_bc(xbd,ybd);
%   input
%          xbd          x boundary coordinate vector
%          ybd          y boundary coordinate vector 
%   IFISS function: DJS; 18 July 2006.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
bc=xbd.*(xbd-1).*ybd.*(ybd-1) + 1;
return