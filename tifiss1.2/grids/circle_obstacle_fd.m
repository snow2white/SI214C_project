function d = circle_obstacle_fd(p)
%CIRCLE_OBSTACLE_FD define circle obstacle geometry for DISTMESH
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao
  g1 = drectangle ( p, 0.0, 8.0,  -1.0, 1.0  );
  g2 = dcircle ( p, 1.75, 0.0, 0.25 );
  d=ddiff(g1,g2);
return
