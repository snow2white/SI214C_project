function d = ell_fd ( p )
%ELL_FD define the L-shape geometry for DISTMESH
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao
  g1 = drectangle ( p, -1.0, 1.0,  0.0, 1.0  );
  g2 = drectangle ( p,  0.0, 1.0, -1.0, 0.0 );
  d = dunion ( g1, g2 );
return