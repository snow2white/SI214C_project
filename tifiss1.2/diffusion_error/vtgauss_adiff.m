function [diffx,diffy] = vtgauss_adiff(s,t,xl,yl)
%VTGAUSS_ADIFF  elementwise permeability field for vector s, t
%   [diffx,diffy] = vtgauss_adiff(s,t,xl,yl);
%   input
%          s         first triangle coordinate
%          t         second triangle coordinate
%          xl        physical element x vertex coordinates 
%          yl        physical element y vertex coordinates  
%
%   calls function: specific_adiff
%    TIFISS function: DJS; 16 January 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao
      nel=length(xl(:,1));
      zero_v = zeros(nel,1); xx=zero_v; yy=xx;
      [xi,dxids,dxidt] = vtshape(s,t);
      for ivtx=1:3
      xx = xx + xi(ivtx) * xl(:,ivtx);
      yy = yy + xi(ivtx) * yl(:,ivtx);
	  end
      [diffx,diffy] = specific_adiff(xx,yy,nel);       
      return