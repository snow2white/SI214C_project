function [a,q,f,ae,qe] = femp1_diff(xy,evt)
%FEMP1_DIFF  set up linear anisotropic diffusion matrices
%   [A,Q,f,Ae,Qe] = femp1_diff(xy,evt);
%   input
%          xy         vertex coordinate vector  
%          evt        element mapping matrix
%   output
%          A          diffusion matrix
%          Q          mass matrix 
%          f          rhs vector
%          Ae         element diffusion matrices
%          Qe         element mass matrices
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%    TIFISS function: DJS; 3 March 2017.
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(evt(:,1));
fprintf('setting up P1 diffusion matrices...  ')
%
% initialise global matrices
      a = sparse(nvtx,nvtx);
      q = sparse(nvtx,nvtx);
      f = zeros(nvtx,1);
%
% inner loop over elements    
        for ivtx = 1:3
        xl_v(:,ivtx) = x(evt(:,ivtx));
        yl_v(:,ivtx) = y(evt(:,ivtx)); 
	    end
        ae = zeros(nel,3,3);
        qe = zeros(nel,3,3);
        fe = zeros(nel,3);  
%
%  One point Gauss rule integration
         sigpt=1/3;
         tigpt=1/3;
         wt=1/2;
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
         rhs = tgauss_source(sigpt,tigpt,xl_v,yl_v);
         [diffx,diffy] = tgauss_adiff(sigpt,tigpt,xl_v,yl_v);        
         for j = 1:3
               for i = 1:3
               ae(:,i,j) = ae(:,i,j) + wt* diffx(:).*dphidx(:,i).*dphidx(:,j) .* invjac(:);
               ae(:,i,j) = ae(:,i,j) + wt* diffy(:).*dphidy(:,i).*dphidy(:,j) .* invjac(:);
               end
           fe(:,j) = fe(:,j) + wt* rhs(:) .* phi(:,j) .* jac(:); %check factor
         end
%
%------------------------------------------------
% 3 point Gauss rule integration
nngpt=3; [s,t,wt]=triangular_gausspoints(nngpt);

%
% loop over Gauss points
      for igpt = 1:nngpt
         sigpt=s(igpt);
         tigpt=t(igpt);
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
         for j = 1:3
               for i = 1:3
               qe(:,i,j) = qe(:,i,j)  + wt(igpt)*phi(:,i).*phi(:,j) .* jac(:);
               end
 	        end
% end of Gauss point loop
      end
%
% perform assembly of global matrix  and source vector
      for krow=1:3
	  nrow=evt(:,krow);	 
          for kcol=1:3
		  ncol=evt(:,kcol);	  
          a = a + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
          q = q + sparse(nrow,ncol,qe(:,krow,kcol),nvtx,nvtx);
          end
%      f(nrow,1) = f(nrow,1) + fe(:,krow)
      for els=1:nel; f(nrow(els),1)=f(nrow(els),1) + fe(els,krow); end
      end
%
fprintf('done\n')
return

