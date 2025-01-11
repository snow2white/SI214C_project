function [bc] = goafem_specific_bc(xbd,ybd)
%GOAFEM_ZERO_BC zero Dirichlet boundary conditions
%  
% [bc] = goafem_specific_bc(xbd,ybd)
%   
% input:
%        xbd    x boundary coordinate vector
%        ybd    y boundary coordinate vector 
%
% output: 
%         bc    boundary conditions on boundary nodes
%
% This function is a copy of original TIFISS function ZERO_BC 
% (DJS; 28 February 2005).
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  bc = zeros(size(xbd));

end % end function