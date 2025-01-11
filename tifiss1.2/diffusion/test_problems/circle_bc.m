function bc = specific_bc(xbd,ybd)
%circle_bc   zero boundary condition on unit circle
%   bc = specific_bc(xbd,ybd);
%   input
%          xbd          x boundary coordinate vector
%          ybd          y boundary coordinate vector 
%   IFISS function: DJS; 28 February 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
bc=0.25*( 1-(xbd.*xbd+ybd.*ybd));
return