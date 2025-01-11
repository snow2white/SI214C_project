function bc = specific_bc(xbd,ybd)
%BC_EX7   nonzero boundary condition for singular solution
%   bc = specific_bc(xbd,ybd);
%   input
%          xbd          x boundary coordinate vector
%          ybd          y boundary coordinate vector 
%   T-IFISS function: AB; 6 October 2017.
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi
  r  = sqrt( xbd.^2 + ybd.^2 );
  theta = atan2(ybd,xbd); 
  bc = (r.^(2/3)) .* sin( (2*theta + pi) ./ 3 );

end % end function
