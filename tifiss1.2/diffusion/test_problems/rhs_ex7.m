function f = specific_rhs(x,y,nel)
%RHS_EX7   RHS forcing function
%   f = specific_rhs(x,y,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements  
%   TIFISS function:  AB; 6 October 2017.
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi
  r = sqrt( x.^2 + y.^2 );
  theta = atan2(y,x); 
  f = (1/3) * (r.^(2/3)) .* sin( (2*theta + pi) ./ 3 );  
end % end function
