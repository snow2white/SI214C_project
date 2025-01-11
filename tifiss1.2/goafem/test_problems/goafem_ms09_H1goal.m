function [g1,g2] = goafem_specific_H1goal(x,y,nel)
%GOAFEM_MS09_H1GOAL Mommer-Stevenson(2009) H1 part of the RHS of the dual problem 
%
% [g1,g2] = goafem_specific_H1goal(x,y,nel)
%   
% input:
%          x    x coordinate vector
%          y    y coordinate vector
%        nel    number of elements
%
% output:
%    [g1,g2]    source vector of the RHS (dual)
% 
% The function computes the source vector \vec{g}=[g1,g2] with
% - g1 = \chi_T,   T=conv{(1,1),(0.5,1),(1,0.5)}
% - g2 = 0
% in the RHS G(v) of the dual problem (see also GOAFEM_FEMP1_ADIFF) 
% as given in [MS09,Example7.3] posed on the square domain (0,1)^2;
%
% Reference:
% [MS09] Mommer, Stevenson, A goal-oriented finite element method with
% convergence rates, SIAM J. Numer. Anal., 47(2)861-866, 2009;
%
% See also GOAFEM_MS09_H1RHS
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  g1 = zeros(nel,1);
  g2 = zeros(nel,1);

  for i = 1:length(x)
      if (0.5 <= x(i) && x(i) <= 1.0 && 1.5-x(i) <= y(i) && y(i) <= 1.0)
          g1(i) = 1.0;
      end
  end

end % end function