function [rhsq,hlsq] = p1res_diff(xy,ev)
%P1RES_DIFF  interior residuals for triangular P1 grid
%   [rhsq,hlsq] = p1res_diff(xy,evt);
%   input
%          xy        vertex coordinate vector  
%          evt       element mapping matrix
%   output
%          rhsq      elementwise L2 residual norms
%          hlsq      elementwise areas
%    PIFISS function: DJS; 31 January 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(ev(:,1));
fprintf('computing P1 interior residuals...  ')
%
%  Mid-side Gauss rule integration
      s(1) = 0.5; t(1) = 0.5;  wt(1)=1/6;
      s(2) = 0.0; t(2) = 0.5;  wt(2)=1/6;
      s(3) = 0.5; t(3) = 0.0;  wt(3)=1/6;
%
         rhsq = zeros(nel,1);
         hlsq = zeros(nel,1);
% inner loop over elements    
         for ivtx = 1:3
         xl_v(:,ivtx) = x(ev(:,ivtx));
         yl_v(:,ivtx) = y(ev(:,ivtx)); 
		 end
% loop over  Gauss points
         for igpt = 1:3
         sigpt=s(igpt);
         tigpt=t(igpt);
         wgt=wt(igpt);
%  evaluate derivatives etc
        [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
         rhs = tgauss_source(sigpt,tigpt,xl_v,yl_v);
         rhsq(:) = rhsq(:) + wgt* rhs(:) .* rhs(:) .* jac(:); 
         hlsq(:) = hlsq(:) + wgt* jac(:);
% end of Gauss point loop
         end
%
fprintf('done\n')
return