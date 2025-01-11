function [diffx,diffy] = specific_adiff(x,y,nel)
%ref_adiff   reference variable diffusion operator 
%   [diffx,diffy] = specific_adiff(x,y,nel);
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements
%   PSFEM function: DJS; 23 June 2006. 
%Copyright (c) 2006 by C.E. Powell and D.J. Silvester (see readme.m)
      diffx =  1*ones(nel,1);  diffy =  1*ones(nel,1); 
      k=find(abs(x)<0.5); diffy(k)=0.01; 
      k=find(abs(y)<0.5); diffx(k)=0.01;  
      return