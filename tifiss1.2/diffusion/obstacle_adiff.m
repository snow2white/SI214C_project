%OBSTACLE_ADIFF   solve anisotropic diffusion problem in obstacle domain
%    TIFISS scriptfile: LR; 24 January 2019.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

%% define geometry
pde=1; domain=17;
obstacle_domain;
load obstacle_grid;
% renumbering nodes for unstructured meshes
[evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);

%% uniform refinement
nref=default('refinement level (default 0)',0);
for i=1:nref
    [xy,evt,bound,eboundt] = p1_refinement(xy,evt,bound,eboundt);
    % renumbering nodes for unstructured meshes
    [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);
    save obstacle_grid.mat evt xy bound eboundt cylinder_choice;
end
% validate mesh
[eex,tve,els] = tedgegen(xy,evt);

%% set up matrices
pmethod=default('P1/P2 approximation 1/2? (default P2)',2);
if pmethod==1
    [A,M,f] = femp1_adiff(xy,evt);  
elseif pmethod==2
    xyp1 = xy;
    boundp1 = bound;
    [xy,evt,bound] = p2grid(xy,evt,bound,eboundt);
    [A,M,f] = femp2_diff(xy,evt,xy,evt);  
end

%% boundary conditions
[Agal,fgal] = nonzerobc(A,f,xy,bound);

%% save resulting system
fprintf('system saved in system_adiff.mat ...\n');
gohome; cd datafiles; save system_adiff pmethod Agal M fgal xy;

%% compute solution
tic;
fprintf('solving linear system ...  ');
x_gal=Agal\fgal;
fprintf('done\n');
etoc=toc; fprintf('Galerkin system solved in %8.3e seconds\n',etoc);
save system_adiff x_gal  -append;

%% plot solution
hplot=1/130;
xmin=min(xy(:,1)); xmax=max(xy(:,1));
ymin=min(xy(:,2)); ymax=max(xy(:,2));
x=xmin:hplot:xmax; y=ymin:hplot:ymax;
tsolplotobs(pmethod,x_gal,evt(:,1:3),xy,x,y,12);
drawnow;

%% compute a posteriori error estimate
if pmethod==1
    pestim=default('error estimation: linear/quadratic bubble functions 1/2 (default 2)',2);  
    errtime=tic;
    if pestim==1
        % linear midpoint hat functions
        fprintf('<strong>diffpost using 3 edge midpoint linear functions...</strong>\n'); 
        [elerr_p,fe,ae] = diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,x_gal,2);
    elseif pestim==2
        % quadratic midpoint bubble functions
        fprintf('\n<strong>diffpost using 4 quadratic bubble functions...</strong>\n');
        [elerr_p,fe,ae] = diffpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,x_gal); 
    else
        error('Value not valid! Choose between either 1 or 2!');
    end
    % boundary correction
    [err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
    fprintf('estimation took %8.3e seconds\n',toc(errtime));
    % plot computed error    
    fprintf('<strong>estimated energy error: %10.4e</strong>\n',norm(elerr_p,2));
    terrplotobsp1(pmethod,x_gal,elerr_p,evt,xy,x,y,101);
elseif pmethod==2
    tic,[elerr_p_p4,fe,ae] = ...
    diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal); toc
    %fprintf('sanity check %10.4e\n',norm(elerr_p_p4-elerr_x_p4,inf))
    terrplotobs(pmethod,x_gal,elerr_p_p4,evt,xy,x,y,102,[xmin,xmax,ymin,ymax]);
end

% end scriptfile
