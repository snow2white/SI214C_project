function [rhs,goal] = goafem_tgauss_L2rhs(s,t,xl,yl)
%GOAFEM_TGAUSS_L2RHS evaluates the L2 part of the RHS of both the primal and dual problems at Gaussian point
%
% [rhs] = goafem_tgauss_L2rhs(s,t,xl,yl)
%
% input:
%                s    reference element x coordinate   
%                t    reference element y coordinate
%               xl    physical element x vertex coordinates 
%               yl    physical element y vertex coordinates
% 
% output:
%       [rhs,goal]    L2 parts of the rhs of the primal and the dual problem
%
% Function(s) called: goafem_specific_L2rhs
%                     goafem_specific_L2goal
%
% See also GOAFEM_TGAUSS_H1RHS
%
%   TIFISS function: MR; 12 June 2018
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
  rhs = goafem_specific_L2rhs(xx,yy,nel);
  
% Dual problem
  goal = goafem_specific_L2goal(xx,yy,nel);
  
end % end function