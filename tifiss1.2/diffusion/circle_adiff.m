%CIRCLE_ADIFF   solve anisotropic diffusion problem in circular domain
%    TIFISS scriptfile: LR; 24 January 2019.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

%% define geometry
pde=1; domain=9;
circle_domain;
load circle_grid;
% renumbering nodes for unstructured meshes
[evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);
% plot the mesh
%figure(10)
imeshplot(xy,evt,bound,eboundt);

% uniform refinement
nref=default('refinement level (default 0)',0);
for i=1:nref    
    [xy,evt,bound,eboundt] = p1_refinement(xy,evt,bound,eboundt);
    [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);
    save circle_grid.mat xy evt bound eboundt;
    %figure(10) 
    %imeshplot(xy,evt,bound,eboundt);
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
xx=xy(:,1); yy=xy(:,2);
xxmin=min(xx); xxmax=max(xx); xxd=(xxmax-xxmin)/256;
x=[xxmin:xxd:xxmax];
yymin=min(yy); yymax=max(yy); yyd=(yymax-yymin)/256;
y=[yymin:yyd:yymax];
tsolplot(pmethod,x_gal,evt(:,1:3),xy,x,y,12,[xxmin,xxmax,yymin,yymax]);
drawnow

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
    terrplot(pmethod,x_gal,elerr_p,evt,xy,x,y,101);
elseif pmethod==2
    tic,[elerr_p_p4,fe,ae] = ...
    diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal); toc
    %fprintf('sanity check %10.4e\n',norm(elerr_p_p4-elerr_x_p4,inf))
    terrplot(pmethod,x_gal,elerr_p_p4,evt,xy,x,y,102,[xxmin,xxmax,yymin,yymax]);
end

% end scriptfile
