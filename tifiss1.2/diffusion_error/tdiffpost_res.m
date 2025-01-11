function [err_res,elerr_res] = tdiffpost_res(jmp,els,rhsq,hlsq)
%TDIFFPOST_RES  computes P1 element residual error estimator
%   [err_res,elerr_res] = tdiffpost_res(jmp,els,rhsq,hlsq);
%   input
%          jmp          elementwise edge flux jumps
%          els          elementwise edge lengths
%          rhsq         elementwise L2 residual norms
%          hlsq         elementwise areas
%   output
%          err_res      global residual error
%          elerr_res    elementwise residual errors
%
%    PIFISS function: DJS; 31 January 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
elerr_res=0.5*sum((els.*els.*jmp.*jmp.*1/4)')' + hlsq.*rhsq;
err_res=sqrt(sum(elerr_res));
elerr_res=sqrt(elerr_res);
fprintf('computing residual error estimator... ')
fprintf('\nestimated global error (in energy):  %10.6e\n',err_res)   
return