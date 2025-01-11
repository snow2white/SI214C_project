function [f1,f2] = goafem_specific_H1rhs(x,y,nel)
%GOAFEM_MOD_MS09_H1RHS Mommer-Stevenson(2009)-type H1 part of the RHS of the primal problem 
%
% [f1,f2] = goafem_specific_H1rhs(x,y,nel)
%   
% input:
%          x    x coordinate vector
%          y    y coordinate vector
%        nel    number of elements
%
% output:
%    [f1,f2]    source vector of the RHS (primal)
% 
% The function computes the source vector \vec{f}=[f1,f2] with
% - f1 = \chi_T,   T=conv{(-1,1),(-1,0.5),(-0.5,1)}
% - f2 = 0
% in the RHS F(v) of the primal problem (see also GOAFEM_FEMP1_ADIFF) 
% as a modification of [MS09,Example7.3].
%
% Reference:
% [MS09] Mommer, Stevenson, A goal-oriented finite element method with
% convergence rates, SIAM J. Numer. Anal., 47(2)861-866, 2009;
%
% See also GOAFEM_MOD_MS09_H1GOAL
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  f1 = zeros(nel,1);
  f2 = zeros(nel,1);

  for i = 1:length(x)
      if (-1 <= x(i) && x(i) <= -0.5 && x(i)+1.5 <= y(i) && y(i) <= 1)
          f1(i) = 1;
      end
  end

end % end function