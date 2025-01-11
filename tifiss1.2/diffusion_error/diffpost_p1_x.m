function [elerr_p,fe,ae] = diffpost_p1_x(xy,evt,eex,tve,els,eboundt,p1sol)
%DIFFPOST_P1_X  legacy code
%   [elerr_p,fe,ae] = diffpost_p1_x(xy,evt,eex,tve,els,eboundt,x_gal);
%   input
%          xy         vertex coordinate vector
%          evt        element mapping matrix
%          eex        presorted element connectivity array
%          tve        presorted edge location array
%          els        element edge size array
%          eboundt    element edge boundary matrix
%          p1sol      P1 solution for diffusion problem
%   output
%          elerr_p    elementwise error estimate
%          fe         elementwise rhs vectors
%          ae         elementwise Poisson problem matrices
%
% calls functions tgauss_source, triangular_gausspoints, gausspoints_oned
%   TIFISS function: DJS; 25 March 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao

fprintf('fast computation of P1 local error estimator...\n')
x=xy(:,1); y=xy(:,2);
nel=length(evt(:,1));
elerr_p=zeros(nel,1);
nnode=4;
%
% set up 7 point Gauss rule with precision 5
nngpt=7; [s,t,wt]=triangular_gausspoints(nngpt); %% nngpt=7;

% inner loop over elements    
    for ivtx = 1:3
        xl_v(:,ivtx) = x(evt(:,ivtx));
        yl_v(:,ivtx) = y(evt(:,ivtx)); 
		end
    ae = zeros(nel,nnode,nnode); elerr=zeros(nnode,nel);
    fe = zeros(nel,nnode);
% loop over Gauss points
    for igpt = 1:nngpt
        sigpt=s(igpt); tigpt=t(igpt); wght=wt(igpt);
% evaluate derivatives etc
        [jac_v,invjac_v,phi_v,dphidx_v,dphidy_v] = tderiv(sigpt,tigpt,xl_v,yl_v);
        [psi_v,dpsidx_v,dpsidy_v] = tqderiv(sigpt,tigpt,xl_v,yl_v);
        rhs_v = tgauss_source(sigpt,tigpt,xl_v,yl_v);
        for j = 1:nnode
             for i = 1:nnode
             ae(:,i,j) = ae(:,i,j)+wght*dpsidx_v(:,i+3).*dpsidx_v(:,j+3).*invjac_v(:);
             ae(:,i,j) = ae(:,i,j)+wght*dpsidy_v(:,i+3).*dpsidy_v(:,j+3).*invjac_v(:);
             end
         fe(:,j) = fe(:,j)  +  wght* rhs_v(:) .* psi_v(:,j+3) .* jac_v(:);
         end
 % end of Gauss point loop
    end
%
%% debug I
res_int=fe; fe = zeros(nel,nnode);

% compute flux jumps
jmp = p1fluxjmps_x(p1sol,eex,tve,xy,evt,eboundt);
% set up one-dimensional Gauss Rule
   nngpt=2; [oneg,onew]=gausspoints_oned(nngpt); %%% nngpt=2
%
% compute residual source vector
        for ig=1:nngpt
            sigpt=oneg(ig); wt=onew(ig);
            for j=1:3
                [s,t] = reorder_s_x(sigpt,evt,j);
                [psi_v,dpsidx_v,dpsidy_v] = vtqderiv(s,t,xl_v,yl_v);
                fe(:,j) = fe(:,j) - wt*(1/2)*jmp(:,j).*psi_v(:,j+3).* els(:,j)./2;
             end
        end
%
% solve for local estimate (sequential code)
%         for ielem = 1:nel
%		  elerr(:,ielem) = squeeze(ae(ielem,1:nnode,1:nnode))\(fe(ielem,1:nnode)');
%         end

%% debug II
res_edge=-fe;
fprintf('internal_res = %7.4e;     edge_res = %7.4e\n',norm(res_int), norm(res_edge));
%disp([res_int])
%disp([res_edge])
fe=res_int-res_edge;
%disp([fe])
                                                    
% quick fix
%AE=ae;         % save element matrices
                                                               
                                                               
% vectorized code
% LDLT factorization
nn=nnode;
dd=zeros(nel,nn); rr=zeros(nel,nn);
for kk=1:nn-1
    for pp=1:kk-1;
        rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,kk,pp);
    end
    dd(1:nel,kk)=ae(1:nel,kk,kk);
    for pp=1:kk-1;
        dd(1:nel,kk)= dd(1:nel,kk)- ae(1:nel,kk,pp).*rr(1:nel,pp);
    end
    for ii=kk+1:nn
        for pp=1:kk-1;
        ae(1:nel,ii,kk)=ae(1:nel,ii,kk)-ae(1:nel,ii,pp).*rr(1:nel,pp);
        end
    ae(1:nel,ii,kk)=ae(1:nel,ii,kk)./dd(1:nel,kk);
    end
end
for pp=1:nn-1;
    rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,nn,pp);
end
dd(1:nel,nn)=ae(1:nel,nn,nn);
for pp=1:nn-1;
    dd(1:nel,nn)= dd(1:nel,nn)- ae(1:nel,nn,pp).*rr(1:nel,pp);
end
% overwrite diagonal entries
for kk=1:nn
    ae(1:nel,kk,kk)= dd(1:nel,kk);
end
%  forward-backward substitutions ...
xx = element_lusolve(ae,fe);
elerr=xx';
                                                               
for ivtx=1:nnode
    elerr_p(:) = elerr_p(:) + fe(:,ivtx) .* elerr(ivtx,:)';
end
elerr_p=sqrt(elerr_p);
fprintf('estimated energy error is %10.4e \n',norm(elerr_p,2))
return
