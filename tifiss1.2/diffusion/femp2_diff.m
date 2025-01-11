function [a,q,f,ae,qe] = femp2_diff(xy,evt,p2xy,p2evt)
%FEMP2_DIFF  set up quadratic anisotropic diffusion matrices
%   [A,Q,f,Ae,Qe] = femp2_diff(xy,evt,p2xy,p2evt);
%   input
%          xy         vertex coordinate vector  
%          evt        element mapping matrix
%          p2xy       P2 node coordinate  vector
%          p2evt      P2 element mapping matrix
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
% Copyright (c) 2011 D.J. Silvester, Qifeng Liao
x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nvtxp2=length(p2xy(:,1));
nel=length(evt(:,1));
fprintf('setting up P2 diffusion matrices...  ')
%
% initialise global matrices
      a = sparse(nvtxp2,nvtxp2);
      q = sparse(nvtxp2,nvtxp2);
      f = zeros(nvtxp2,1);
%
% inner loop over elements    
        for ivtx = 1:3
        xl_v(:,ivtx) = x(evt(:,ivtx));
        yl_v(:,ivtx) = y(evt(:,ivtx)); 
	    end
        ae = zeros(nel,6,6);
        qe = zeros(nel,6,6);
        fe = zeros(nel,6);  

%------------------------------------------------
% 7 point Gauss rule integration
nngpt=7;
[s,t,wt]=triangular_gausspoints(nngpt);
            
for igpt=1:nngpt
       sigpt=s(igpt);
       tigpt=t(igpt);
       wtigpt=wt(igpt);
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
         [psi,dpsidx,dpsidy] = tqderiv(sigpt,tigpt,xl_v,yl_v);
         rhs = tgauss_source(sigpt,tigpt,xl_v,yl_v);
         [diffx,diffy] = tgauss_adiff(sigpt,tigpt,xl_v,yl_v);        
         for j = 1:6
               for i = 1:6
               ae(:,i,j) = ae(:,i,j) + wtigpt* diffx(:).*dpsidx(:,i).*dpsidx(:,j) .* invjac(:);
               ae(:,i,j) = ae(:,i,j) + wtigpt* diffy(:).*dpsidy(:,i).*dpsidy(:,j) .* invjac(:);
               end
               fe(:,j) = fe(:,j) + wtigpt* rhs(:) .* psi(:,j) .* jac(:); %check factor
         end
         for j = 1:6
               for i = 1:6
               qe(:,i,j) = qe(:,i,j)  + wtigpt*psi(:,i).*psi(:,j) .* jac(:);
               end
 	        end
% end of Gauss point loop
end
%
% perform assembly of global matrix  and source vector
      for krow=1:6
	     nrow=p2evt(:,krow);	 
          for kcol=1:6
		  ncol=p2evt(:,kcol);	  
          a = a + sparse(nrow,ncol,ae(:,krow,kcol),nvtxp2,nvtxp2);
          q = q + sparse(nrow,ncol,qe(:,krow,kcol),nvtxp2,nvtxp2);
          end
          for els=1:nel; f(nrow(els),1)=f(nrow(els),1) + fe(els,krow); end
      end
%
fprintf('done\n')
return

