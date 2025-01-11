%UNITCRACK outlines unit crack domain (0,1)^2 \(1/2,1)x{1/2} (crack on the right)
%   TIFISS scriptfile: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi
 
  global slope1crack
  global slope2crack
  hold on  
  
% Lower diagonal line
  xP1 = 1; 
  yP1 = 0.5 + slope1crack*(xP1-0.5);
  plot([0.5,xP1],[0.5,yP1],'-k');
% Upper diagonal line
  xP2 = 1; 
  yP2 = 0.5 + slope2crack*(xP2-0.5);
  plot([0.5,xP2],[0.5,yP2],'-k');
% Rest of lines       
  plot([0,1],   [0,0], '-k');
  plot([1,xP1], [0,yP1], '-k');
  plot([xP2,1], [yP2,1], '-k');
  plot([1,0], [1,1], '-k');
  plot([0,0], [1,0], '-k');
  axis square;
  hold off;
  
% end scriptfile