function [diffx,diffy] = specific_adiff(x,y,nel)
%analytic_adiff   analytic unit circle diffusion operator 
%   [diffx,diffy] = specific_adiff(x,y,nel);
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements
%   PSFEM function: DJS; 19 July 2006. 
%Copyright (c) 2006 by C.E. Powell and D.J. Silvester (see readme.m)
      diffx =  1*ones(nel,1);  diffy =  1*ones(nel,1); 
      r=sqrt(x.*x+y.*y);
      k=find(r<0.7); diffx(k)=100;  diffy(k)=100; 
      k=find(r<0.3); diffx(k)=1/10; diffy(k)=1/10;   
      return