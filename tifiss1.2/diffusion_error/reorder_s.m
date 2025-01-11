function [s,t] = reorder_s(sigpt,edge_number)
%REORDER_S locate 1D Gauss points on the reference element edges
% [s,t] = reorder_s(sigpt,edge_number)
%   input
%         sigpt        1D Gauss Point
%         edge_number  edge number
%   output
%         s            x-coordinate of Gauss point on reference element
%         t            y-coordinate of Gauss point on the reference element
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

%s = (1 + sigpt)/2.0;   % Map form [-1,1] to [0,1]

  if     edge_number == 1
         % first DIAGONAL reference edge
         s = (1.0 + sigpt)/2.0;  
         t = 1 - s;
  elseif edge_number == 2
         % second LEFT reference edge
         s = 0.0;
         t = (1.0 + sigpt)/2.0;
  elseif edge_number == 3
         % third BOTTOM reference edge
         s = (1.0 + sigpt)/2.0;  
         t = 0.0;
  end


end  % end function

