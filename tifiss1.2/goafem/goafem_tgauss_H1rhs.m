function [f1,f2,g1,g2] = goafem_tgauss_H1rhs(s,t,xl,yl)
%GOAFEM_TGAUSS_H1RHS evaluates the H1 part of the RHS of both the primal and dual problems at Gaussian point
%
% [f1,f2,g1,g2] = goafem_tgauss_H1rhs(s,t,xl,yl)
%
% input:
%                s    reference element x coordinate   
%                t    reference element y coordinate
%               xl    physical element x vertex coordinates 
%               yl    physical element y vertex coordinates
% 
% output:
%          [f1,f2]    H1 part of the rhs of the primal problem
%          [g1,g2]    H1 part of the rhs of the dual problem
%
% Function(s) called: goafem_specific_H1rhs
%                     goafem_specific_H1goal
%
% See also GOAFEM_TGAUSS_L2RHS 
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel      = length(xl(:,1));
  xx       = zeros(nel,1);
  yy       = xx;
  [xi,~,~] = tshape(s,t);
  
  for ivtx = 1:3 
      xx = xx + xi(ivtx) * xl(:,ivtx);
      yy = yy + xi(ivtx) * yl(:,ivtx);
  end
  
% Primal problem
  [f1,f2] = goafem_specific_H1rhs(xx,yy,nel);

% Dual problem
  [g1,g2] = goafem_specific_H1goal(xx,yy,nel);

end % end function