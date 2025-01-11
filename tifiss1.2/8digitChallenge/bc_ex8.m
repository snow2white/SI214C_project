function bc = specific_bc(xbd,ybd)
%BC_EX8   nonzero boundary condition for singular solution
%   bc = specific_bc(xbd,ybd);
%   input
%          xbd          x boundary coordinate vector
%          ybd          y boundary coordinate vector 
%   T-IFISS function: AB; 6 October 2017.
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  bc = (1-xbd).^2;

end % end function
