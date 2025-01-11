function [diffx,diffy] = specific_adiff(x,y,nel)
%ADIFF_EX7 variable diffusion operator 
%   [diffx,diffy] = specific_adiff(x,y,nel);
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements
%   output
%          diffx      first component 
%          diffy      second component
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

   diffx = 1 - (x.^2 + y.^2)/4; 
   diffy = 1 - (x.^2 + y.^2)/4;

end  % end function
