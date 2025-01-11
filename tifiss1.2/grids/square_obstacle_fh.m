function d = square_obstacle_fh(p)
%SQUARE_OBSTACLE_FH define square obstacle geometry for DISTMESH
%    TIFISS function: DJS; 12 February 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao
  d = 0.2+0.3* drectangle ( p, 1.75,2.25, -0.25, 0.25 );
return
