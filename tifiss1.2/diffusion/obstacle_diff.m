%OBSTACLE_DIFF   legacy code
%    TIFISS scriptfile: DJS; 5 March 2017.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

%% define geometry
pde=1; domain=17;
obstacle_domain
load obstacle_grid

% uniform refinement
nref=default('refinement level (default 0)',0);
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
    xyp1=xy;
    [xy,evt,bound]=p2grid(xy,evt,bound,eboundt);
    [A,M,f] = femp2_diff(xy,evt,xy,evt);
end 
%% boundary conditions
   [Agal,fgal] = nonzerobc(A,f,xy,bound);
%% save resulting system
fprintf('system saved in system_adiff.mat ...\n')
gohome
cd datafiles
save system_adiff qmethod Agal M  fgal  xy
%% compute solution
tic
fprintf('solving linear system ...  ')
x_gal=Agal\fgal;
fprintf('done\n')
etoc=toc; fprintf('Galerkin system solved in  %8.3e seconds\n',etoc) 
save system_adiff x_gal  -append
%
%% plot solution
hplot=1/300;
xmin=min(xy(:,1)); xmax=max(xy(:,1));
ymin=min(xy(:,2)); ymax=max(xy(:,2));
x=xmin:1/300:xmax;
y=ymin:1/300:ymax;
tsolplotobs(qmethod,x_gal,evt(:,1:3),xy,x,y,12);
drawnow

%% compute a posteriori error estimate
if qmethod==1
   tic,[elerr_p,fe,ae] = diffpost_p1_x(xy,evt,eex,tve,els,eboundt,x_gal);toc
   [elerror_p] = diffpost_p1_bc_x(ae,fe,xy,evt,eboundt);
   terrplotobs(qmethod,x_gal,elerror_p,evt,xy,x,y,101,[xmin,xmax,ymin,ymax]);
   fprintf('<strong>Estimated energy error: %10.4e</strong>\n',norm(elerror_p,2));
elseif qmethod==2
   tic,[elerr_p_p4,fe,ae] = ...
   diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal); toc
   %fprintf('sanity check %10.4e\n',norm(elerr_p_p4-elerr_x_p4,inf))
   terrplotobs(qmethod,x_gal,elerr_p_p4,evt,xy,x,y,102,[xmin,xmax,ymin,ymax]);
end
