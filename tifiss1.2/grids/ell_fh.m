function h = ell_fh ( p )
%ELL_FH set mesh density of L-shape geometry for DISTMESH
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao

global mesh_density_ratio
  x=p(:,1); y=p(:,2);
  n=length(x);
  h = (x.^2+y.^2)+mesh_density_ratio;  
return