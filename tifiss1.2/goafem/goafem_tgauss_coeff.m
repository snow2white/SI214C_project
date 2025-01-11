function [diffx,diffy] = goafem_tgauss_coeff(s,t,xl,yl)
%GOAFEM_TGAUSS_COEFF evaluates diffusion coefficient at triangle Gauss point
%
% [diffx,diffy] = goafem_tgauss_coeff(s,t,xl,yl)
%
% input:
%              s    reference element x coordinate   
%              t    reference element y coordinate
%             xl    physical element x vertex coordinates 
%             yl    physical element y vertex coordinates
%
% output
%  [diffx,diffy]    diffusion coefficient at Gauss point
% 
%
% This is a copy of original PIFISS function TGAUSS_ADIFF (DJS; 31 January 2007)
%
% Function(s) called: goafem_specific_coeff
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel      = length(xl(:,1));
  xx       = zeros(nel,1); 
  yy       = xx;
  [xi,~,~] = tshape(s,t);
  
  for ivtx=1:3
      xx = xx + xi(ivtx) * xl(:,ivtx);
      yy = yy + xi(ivtx) * yl(:,ivtx);
  end
  
  [diffx,diffy] = goafem_specific_coeff(xx,yy,nel);

end % end function