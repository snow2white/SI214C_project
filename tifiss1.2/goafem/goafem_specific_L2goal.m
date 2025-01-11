function [g] = goafem_specific_L2goal(x,y,nel)
%GOAFEM_PO09_L2GOAL Prudhomme-Oden mollifier L2 part of RHS of the dual problem
%   
% [g] = goafem_specific_L2goal(x,y,nel)
%
% input:
%        x    x coordinate vector
%        y    y coordinate vector
%      nel    number of elements
%
% outpu:
%        g    Prudhomme-Oden mollifier 
%
% The function computes the deterministic L2 part g0(x) of RHS G(v) of the dual 
% problem (see also GOAFEM_FEMP1_ADIFF) defined as the "mollifier" as 
% given in:
% [PO99] Prudhomme, Oden, On goal-oriented error estimation for elliptic 
% problems: application to the control of pointwise errors, Comput.
% Methods. Appl. Mech. Engrg. 176(1-4):313-331, 1999.
%
% See also GOAFEM_MOLLIFIER
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri
   
% Load constant parameters for pointwise estimation 
  load constmollifier_po99.mat;
   
  g = zeros(nel,1);
   
  for i = 1:length(x)
      % Support of the Ball B((x0,y0);r)
      dp = sqrt( (x0 - x(i))^2 + (y0 - y(i))^2 );
      if dp < r
          g(i) = C * exp( -r^2 / (r^2 - dp^2) ); 
      end
  end 
  
end % end function