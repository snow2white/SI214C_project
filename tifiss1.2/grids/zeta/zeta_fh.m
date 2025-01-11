function [h] = zeta_fh(p)
%ZETA_FH sets mesh density of the Z-shaped geometry for DISTMESH
%
%   TIFISS function: AB; 07 June 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi


  global mesh_density_ratio
  
  x = p(:,1);
  y = p(:,2);
  
  % Circular density
  h = (x.^2 + y.^2) + mesh_density_ratio; 
  
  % Uniform densitity: uncomment the circular density and use uniform one 
%   h = ones(size(p,1),1) + mesh_density_ratio;
  
end % end function