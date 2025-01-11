function ff = tgauss_source(s,t,xl,yl)
%TGAUSS_SOURCE  evaluates source term at triangle Gauss point
%   ff = tgauss_source(s,t,xl,yl);
%   input
%          s         reference element x coordinate   
%          t         reference element y coordinate
%          xl        physical element x vertex coordinates 
%          yl        physical element y vertex coordinates  
%
%   calls function: specific_rhs
%    PIFISS function: DJS; 31 January 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
      nel=length(xl(:,1));
      zero_v = zeros(nel,1); xx=zero_v; yy=xx;
      [xi,dxids,dxidt] = tshape(s,t);
      for ivtx=1:3 
      xx = xx + xi(ivtx) * xl(:,ivtx);
      yy = yy + xi(ivtx) * yl(:,ivtx);
	  end
      ff=specific_rhs(xx,yy,nel);
      return