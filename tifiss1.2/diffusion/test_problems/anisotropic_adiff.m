function [diffx,diffy] = specific_adiff(x,y,nel)
%anisotropic_adiff   anisotropic diffusion operator
%   [diffx,diffy] = specific_adiff(x,y,nel);
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements
%   PSFEM function: DJS; 23 June 2006. 
%Copyright (c) 2006 by C.E. Powell and D.J. Silvester (see readme.m)
      diffx =  1*ones(nel,1);  diffy =  100*ones(nel,1); 
      return
