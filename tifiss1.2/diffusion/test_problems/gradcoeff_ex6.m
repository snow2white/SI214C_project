function [dcoeffdx,dcoeffdy] = specific_gradcoeff(x,y,nel)
%GRADCOEFF_EX6 Computes the gradient of the coefficients a(x,y)
%
%   [dcoeffdx,dcoeffdy] = gradcoeff_ex6(x,y,nel)
%
%   input:
%                  x     x coordinate vector
%                  y     y coordinate vector 
%                nel     number of elements
%   output:
%           dcoeffdx     first component derivative
%           dcoeffdy     second component derivative
%
% SEE ALSO ADIFF_EX6
%
% Last update 07/03/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi


% NOTE: 
% -------------------------------------------------------------------------
% Choose either the 1st, 2nd, or 3rd gradient of the coefficient
% accordingly to what is selected in  
% diffusion/test_problems/adiff_ex6.m function


% Gradient first example (default one):
  lambda = 2.183365648442389;
  omega_x = 0.653271187094403;      alpha_x = 0.75835723479235;
  omega_y = 0.653271187094403;      alpha_y = 0.75835723479235;
  type_x = -3;                      type_y = 2;

  fx = alpha_x * ( ((1 + type_x)/2) * cos(x*omega_x) + ((1 - type_x)/2) * sin(x*omega_x) );
  fy = alpha_y * ( ((1 + type_y)/2) * cos(y*omega_y) + ((1 - type_y)/2) * sin(y*omega_y) );
       
  dfxdx = alpha_x * omega_x * ( -((1 + type_x)/2)*sin(x*omega_x) + ((1 - type_x)/2)*cos(x*omega_x) );
  dfydy = alpha_y * omega_y * ( -((1 + type_y)/2)*sin(y*omega_y) + ((1 - type_y)/2)*cos(y*omega_y) ); 
        
  dcoeffdx = 0.2 * sqrt(3.0e0) * sqrt(lambda) * dfxdx .* fy;  
  dcoeffdy = 0.2 * sqrt(3.0e0) * sqrt(lambda) * fx .* dfydy;  
      
   
  
% % Gradient second example:
%   dcoeffdx = - (0.2) * (0.65) * ( sin(0.65*x) .* cos(0.65*y) );  
%   dcoeffdy = - (0.2) * (0.65) * ( cos(0.65*x) .* sin(0.65*y) );
      
      
    
% % Gradient third example:
%   dcoeffdx = - x;  
%   dcoeffdy = - y;


end   % end function