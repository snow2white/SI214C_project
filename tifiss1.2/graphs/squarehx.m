function squarehx(rd_type)
%SQUAREHX outlines the square domain with central square hole
%
% input:
%          rd_type     example random coefficient
%
% According to rd_type, either the domain in [0,1]^2 or in [-1,1]^2 is plotted
%
%   TIFISS function: AB; 05 June 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi
 
  hold on;

  if rd_type == 1
      pOne = 2 * 0.375 - 1; 
      pTwo = 2 * 0.625 - 1;
      % External square
      plot([-1,1],[-1,-1],'-k');
      plot([1,1],[-1,1],'-k');
      plot([1,-1],[1,1],'-k');
      plot([-1,-1],[1,-1],'-k');
      % Internal square  
      plot([pOne,pTwo],[pOne,pOne],'-k');
      plot([pTwo,pTwo],[pOne,pTwo],'-k');
      plot([pTwo,pOne],[pTwo,pTwo],'-k');
      plot([pOne,pOne],[pTwo,pOne],'-k');
  else
      % External square
      plot([0,0],[1,0],'-k');
      plot([1,0],[1,1],'-k');
      plot([1,1],[0,1],'-k');
      plot([0,1],[0,0],'-k');
      % Internal square  
      plot([0.375,0.625],[0.375,0.375],'-k');
      plot([0.625,0.625],[0.375,0.625],'-k');
      plot([0.625,0.375],[0.625,0.625],'-k');
      plot([0.375,0.375],[0.625,0.375],'-k');
  end
  
  axis square;

  hold off;

end % end function