function d = circle_obstacle_fh(p)
%CIRCLE_OBSTACLE_FH define circle obstacle density for DISTMESH
%    TIFISS function: DJS; 12 February 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao
d = 0.1+0.3*dcircle ( p, 1.75, 0.0, 0.25 );
return
