function d = square_obstacle_fd(p)
%SQUARE_OBSTACLE_FD define square obstacle geometry for DISTMESH
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao
  g1 = drectangle ( p, 0.0, 8.0,  -1.0, 1.0  );
  g2 = drectangle ( p, 1.75,2.25, -0.25, 0.25 );
  d=ddiff(g1,g2);
return
