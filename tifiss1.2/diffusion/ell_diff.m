%ELL_DIFF  legacy code
%   TIFISS scriptfile: DJS; 5 March 2017.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao
clear 
%% define geometry
pde=1; domain=2;
mesh_type=default('structured/unstructured mesh : 1/2 (default 2)',2);
if mesh_type==1
ell_domain
load ell_grid
[evt,eboundt] = p1grid(xy,mv,bound,mbound);
elseif mesh_type==2
ell_domain_unstructured
load ell_grid
figure(10)
imeshplot(xy,evt,bound,eboundt);
else
    error('Oops. invalid mesh type')
end

nref=default('refinement level? (default 0)',0);
for i=1:nref,
    [xy,evt,bound,eboundt] = p1_refinement(xy,evt,bound,eboundt);
    end
% validate mesh
[eex,tve,els] = tedgegen(xy,evt);

%% set up matrices
qmethod=default('P1/P2 approximation 1/2? (default P2)',2);
if qmethod==1
   [A,M,f] = femp1_diff(xy,evt);
elseif qmethod==2
    [xy,evt,bound]=p2grid(xy,evt,bound,eboundt);
    [A,M,f] = femp2_diff(xy,evt,xy,evt);
end 
%% boundary conditions
[Agal,fgal] = nonzerobc(A,f,xy,bound);
%
%% save resulting system
fprintf('system saved in system_adiff.mat ...\n')
gohome
cd datafiles
save system_adiff qmethod Agal M  fgal  xy
%% compute solution
tic, fprintf('solving linear system ...  ')
x_gal=Agal\fgal;
fprintf('done\n')
etoc=toc; fprintf('Galerkin system solved in  %8.3e seconds\n',etoc)
save system_adiff x_gal  -append

%% plot solution
xx=xy(:,1); yy=xy(:,2);
xxmin=min(xx); xxmax=max(xx); xxd=(xxmax-xxmin)/256;
x=[xxmin:xxd:xxmax];
yymin=min(yy); yymax=max(yy); yyd=(yymax-yymin)/256;
y=[yymin:yyd:yymax];
tsolplot(qmethod,x_gal,evt,xy,x,y,12,domain);
drawnow
%
%% compute a posteriori error estimate
if qmethod==1
   tic,[elerr_p,fe,ae] = diffpost_p1_x(xy,evt,eex,tve,els,eboundt,x_gal);toc
   [elerror_p] = diffpost_p1_bc_x(ae,fe,xy,evt,eboundt);
   terrplotl(qmethod,x_gal,elerror_p,evt,xy,x,y,101);
   fprintf('<strong>Estimated energy error: %10.4e</strong>\n',norm(elerror_p,2));
elseif qmethod==2
   tic,[elerr_p_p4,fe,ae] = ...
        diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal); toc
   %fprintf('sanity check %10.4e\n',norm(elerr_p_p4-elerr_x_p4,inf))
   terrplotl(qmethod,x_gal,elerr_p_p4,evt,xy,x,y,102);
end
    
