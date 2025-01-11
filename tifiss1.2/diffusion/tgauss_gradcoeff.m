function [dcoeffdx,dcoeffdy] = tgauss_gradcoeff(s,t,xl,yl)
%TGAUSS_GRADCOEFF evaluates permeability field at triangle Gauss point
%   [dcoeffdx,dcoeffdy] = tgauss_gradcoeff(s,t,xl,yl)
%   input
%          s          reference element x coordinate   
%          t          reference element y coordinate
%          xl         physical element x vertex coordinates 
%          yl         physical element y vertex coordinates  
%   output
%          dcoeffdx   first component derivative
%          dcoeffdy   second component derivative
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  nel = length(xl(:,1));
  zero_v = zeros(nel,1); 
  xx = zero_v;
  yy = xx;
  
  [phi_e,dphids,dphidt] = tshape(s,t);
  
  for ivtx=1:3
      xx = xx + phi_e(ivtx) * xl(:,ivtx);
      yy = yy + phi_e(ivtx) * yl(:,ivtx);
  end  
  
  [dcoeffdx,dcoeffdy] = specific_gradcoeff(xx,yy,nel);

end   % end function
