function [dcoeffdx,dcoeffdy] = goafem_tgauss_gradcoeff(s,t,xl,yl)
%GOAFEM_TGAUSS_GRADCOEFF evaluates derivatives of the coefficient at triangle Gauss point
%   
% [dcoeffdx,dcoeffdy] = goafem_tgauss_gradcoeff(s,t,xl,yl)
%
% input:
%           s   reference element x coordinate   
%           t   reference element y coordinate
%          xl   physical element x vertex coordinates 
%          yl   physical element y vertex coordinates  
% output:
%    dcoeffdx   first component derivative
%    dcoeffdy   second component derivative
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nel         = length(xl(:,1));
  xx          = zeros(nel,1);
  yy          = xx;
  [phi_e,~,~] = tshape(s,t);
  
  for ivtx=1:3
      xx = xx + phi_e(ivtx) * xl(:,ivtx);
      yy = yy + phi_e(ivtx) * yl(:,ivtx);
  end  
  
  [dcoeffdx,dcoeffdy] = goafem_specific_gradcoeff(xx,yy,nel);

end % end function
