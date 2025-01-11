function [elerror_p] = diffpost_p1_bc_x(aez,fez,xy,evt,eboundt)
%DIFFPOST_P1_BC_X legacy code
% [error_p,elerror_p,aez] = diffpost_p1_bc_x(aez,fez,xy,evt,eboundt);
%   input
%          aez       elementwise Poisson problem matrices (LDLT factored)
%          fez       elementwise rhs vectors
%          xy        vertex coordinate vector  
%          ev        element mapping matrix
%          eboundt   element edge boundary matrix 
%   output
%          elerror_p   elementwise error estimate
%   calls function: localbc_p
%   TIFISS function: DJS; 25 March 2016
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao
      x=xy(:,1); y=xy(:,2);
      nel=length(evt(:,1));
      lev=[evt,evt(:,1)]; elerror_p=zeros(nel,1);
%
% recompute contributions from elements with Dirichlet boundaries
      nbde=length(eboundt(:,1));
      ebdy = zeros(nel,1);
      edge = zeros(nel,1);
% isolate boundary elements
      for el = 1:nbde
      ee = eboundt(el,1);
      ebdy(ee) = ebdy(ee)+1; edge(ee)=eboundt(el,2);
      end  
%
% two edge elements
      k2=find(ebdy==2);
      nel2b=length(k2);
% loop over two edge elements
      for el = 1:nel2b
      el2e=k2(el);
      kk=find(eboundt(:,1) == el2e);
      edges=eboundt(kk,2);
% set up original matrix and RHS vector
	  ae=squeeze(aez(el2e,1:4,1:4)); 
      fe=fez(el2e,:)';
% reconstruct original element matrix
      lae=eye(4,4)+tril(ae,-1); dae=diag(diag(ae));
      ae=lae*dae*lae';
% set up local coordinates and impose interpolated error as Dirichlet bc
      xl=x(lev(el2e,:)); yl=y(lev(el2e,:)); 
      [bae,fe] = localbc_p(ae,fe,edges,xl,yl);
     fprintf('\n<strong>Warning:</strong> element %g has two boundary edges\n',el2e)
% refactorise the local problem and store
       [lae,dae]=ldl(bae);
       bae = lae+lae'+dae-2*eye(4,4);
       aez(el2e,1:4,1:4)= bae; fez(el2e,:)=fe;
      end
% end of element loop
%
% one edge elements
      k1=find(ebdy==1);
      nel1b=length(k1);
% loop over one edge elements
      for el = 1:nel1b
      el1e=k1(el);
      kk=find(eboundt(:,1) == el1e);
      edges=eboundt(kk,2);
% set up original matrix and RHS vector 
      fe=fez(el1e,:)';
	  ae=squeeze(aez(el1e,1:4,1:4));
% reconstruct original element matrix
      lae=eye(4,4)+tril(ae,-1); dae=diag(diag(ae));
      ae=lae*dae*lae';
% set up local coordinates and impose interpolated error as Dirichlet bc
      xl=x(lev(el1e,:)); yl=y(lev(el1e,:));
      [bae,fe] = localbc_p(ae,fe,edges,xl,yl);
% refactorise the local problem and store
      [lae,dae]=ldl(bae);
      bae = lae+lae'+dae-2*eye(4,4);
      aez(el1e,1:4,1:4)= bae; fez(el1e,:)=fe;
      end
% end of element loop
%
%  forward-backward substitutions ...
      xx = element_lusolve(aez,fez);
      elerr=xx';
        
      for ivtx=1:4
      elerror_p(:) = elerror_p(:) + fez(:,ivtx) .* elerr(ivtx,:)';
      end
      elerror_p = sqrt(elerror_p);
      fprintf('boundary correction done\n')
 return
