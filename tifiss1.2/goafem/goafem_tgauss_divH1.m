function [divf,divg] = goafem_tgauss_divH1(s,t,xl,yl)
%GOAFEM_TGAUSS_DIVH1 evaluates the divergence of the H1 part of the RHS of both primal and dual problems at Gaussian point
%
% [divf,divg] = goafem_tgauss_divH1(s,t,xl,yl)
%
% input: 
%                s    reference element x coordinate   
%                t    reference element y coordinate
%               xl    physical element x vertex coordinates 
%               yl    physical element y vertex coordinates
% 
% output:
%      [divf,divg]    divergence of source vector of RHS (primal and dual)
%
% Function(s) called: goafem_specific_divH1rhs
%                     goafem_specific_divH1goal
%
% See also GOAFEM_TGAUSS_H1RHS
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
  divf = goafem_specific_divH1rhs(xx,yy,nel);

% Dual problem
  divg = goafem_specific_divH1goal(xx,yy,nel);
  
end % end function