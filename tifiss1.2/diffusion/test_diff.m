%TEST_DIFF   solve isotropic problem in half L-shaped domain
%   TIFISS scriptfile: DJS; 21 December 2017.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

%% define geometry
pde=1; domain=12;
test_domain
load test_grid
figure(10)
imeshplot(xy,evt,bound,eboundt);
% validate mesh
[eex,tve,els] = tedgegen(xy,evt);
pause

%------------------- uniform refinement
nref=default('refinement level (default 0)',0);
for i=1:nref,
[xy,evt,bound,eboundt] = p1_refinement(xy,evt,bound,eboundt);
end
[eex,tve,els] = tedgegen(xy,evt); % validate mesh

qmethod=default('P1/P2 approximation 1/2? (default P1)',1);
if qmethod==1,
%-------------------P1 approximation
[A,M,f,Ae,Qe] = femp1_diff(xy,evt);
[Agal,fgal] = nonzerobc(A,f,xy,bound);
%
tic
fprintf('solving linear system ...  ')
x_gal=Agal\fgal;
fprintf('done\n')
etoc=toc; fprintf('Galerkin system solved in  %8.3e seconds\n',etoc)
% plot solution
tsolplotx(x_gal,xy,evt,12)
drawnow,
%tic,[elerr_x,fe,ae] = diffpost_p1_x(xy,evt,eex,tve,els,eboundt,x_gal); toc
tic; [elerr_p,fe,ae] = diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,x_gal,1); toc;
%tic; [elerr_p,fe,ae] = diffpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,x_gal); toc;
%fprintf('sanity check %10.4e\n',norm(elerr_p-elerr_x,inf))
[err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
fprintf('<strong>Estimated energy error: %10.4e</strong>\n',norm(elerr_p,2));
elseif qmethod==2
%-------------------P2 approximation
[xy,evt,bound]=p2grid(xy,evt,bound,eboundt);
[A,M,f,Ae,Qe] = femp2_diff(xy,evt,xy,evt);
[Agal,fgal] = nonzerobc(A,f,xy,bound);
%
tic
fprintf('solving linear system ...  ')
x_gal=Agal\fgal;
fprintf('done\n')
etoc=toc; fprintf('Galerkin system solved in  %8.3e seconds\n',etoc)
% plot solution
tsolplotx(x_gal,xy,evt,22)
drawnow
tic,[elerr_p_p4,fe,ae] = ...
diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal); toc
fprintf('<strong>Estimated energy error: %10.4e</strong>\n',norm(elerr_p_p4,2));
else
error('Oops.. approximation undefined!')
end
