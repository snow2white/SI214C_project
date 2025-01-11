function [diffx,diffy] = specific_adiff(x,y,nel)
%ADIFF_EX6 Variable diffusion operator 
%
%   [diffx,diffy] = adiff_ex6(x,y,nel)
%
%   input:
%               x      x coordinate vector
%               y      y coordinate vector 
%             nel      number of elements
%   output:
%           diffx      first component 
%           diffy      second component
%
% SEE ALSO GRADCOEFF_EX6
% 
% Last update 07/03/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

% NOTE: 
% -------------------------------------------------------------------------
% Choose either the 1st, 2nd, or 3rd coefficient example and 
% uncomment accordingly the corresponding gradient in 
% diffusion/test_problems/gradcoeff_ex6.m function

  
% First example (default one):
  lambda = 2.183365648442389;
  omega_x = 0.653271187094403;    alpha_x = 0.75835723479235;
  omega_y = 0.653271187094403;    alpha_y = 0.75835723479235;
  type_x = -3;                    type_y = 2;
 
  fx = alpha_x * ( ((1 + type_x)/2) * cos(x*omega_x) + ((1 - type_x)/2) * sin(x*omega_x) );
  fy = alpha_y * ( ((1 + type_y)/2) * cos(y*omega_y) + ((1 - type_y)/2) * sin(y*omega_y) );
   
  a = ones(nel,1);  
  a(:,1) = a(:,1) + 0.2 * sqrt(3.0e0) * sqrt(lambda) * fx .* fy;
  
   
  
% % Second example:
%   a = ones(nel,1);  
%   a(:,1) = a(:,1) + 0.2 * cos(0.65*x) .* cos(0.65*y);
  
  
% % Third example:
%   a = ones(nel,1);  
%   a(:,1) = a(:,1) - 0.5 * ((x .* x) + (y .* y));
  
  
% Assign coefficient  
  diffx = a(:,1);
  diffy = a(:,1);

end  % end function