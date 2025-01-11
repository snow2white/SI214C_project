%ZETAX outlines Z-shaped domain
%
%   TIFISS scriptfile: AB; 07 June 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi
 
  global slope
  
  hold on
  
  % Diagonal edge
  if slope > 1
      xp = -1/slope;  
      yp = -1;
      plot([0,xp],[0,yp],'-k');
      plot([xp 0.0],[yp -1.0],'-k');     
  elseif slope == 1
      xp = -1.0;    
      yp = -1.0;
      plot([0,xp],[0,yp],'-k');
      plot([xp 0.0],[yp -1.0],'-k');
  else%if slope < 1
      xp = -1;        
      yp = -slope; 
      plot([0,xp],[0,yp],'-k');
      plot([xp -1.0],[yp -1.0],'-k')
      plot([-1.0 0.0],[-1.0 -1.0],'-k');
  end
  
  plot([0,1],[-1,-1],'-k');
  plot([1,1],[-1,1],'-k');
  plot([1,-1],[1,1],'-k');
  plot([-1,-1],[1,0],'-k');
  plot([-1,0],[0,0],'-k');

  hold off

% end script