%ELL_ADIFF  solve variable diffusion problem in L-shaped domain
%   TIFISS scriptfile: LR; 24 January 2019.
% Copyright (c) 2017 D.J. Silvester, Alex Bespalov, Qifeng Liao

clear 
%% define geometry
pde=1;
domain=2;
mesh_type = default('structured/unstructured mesh : 1/2 (default 1)',1);

if mesh_type==1
    ell_domain
    load ell_grid
    [evt,eboundt] = p1grid(xy,mv,bound,mbound);
elseif mesh_type==2
    ell_domain_unstructured
    load ell_grid
    % renumbering nodes for unstructured meshes
    [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);
    %figure(10)
    imeshplot(xy,evt,bound,eboundt);
else
    error('Oops. invalid mesh type')
end

% Choose red or bisec3 uniform sub-division for spatial estimation 1/2
subdivPar = 2;

nref=default('refinement level? (default 0)',0);
for i=1:nref,
    [xy,evt,bound,eboundt] = p1_refinement(xy,evt,bound,eboundt);
    [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);
    save circle_grid.mat xy evt bound eboundt;
    %figure(10) 
    %imeshplot(xy,evt,bound,eboundt);
end

% validate mesh
[eex,tve,els] = tedgegen(xy,evt);

%% set up matrices
pmethod = default('P1/P2 approximation 1/2? (default P1)',1);
if pmethod==1
   [A,M,f] = femp1_adiff(xy,evt);
elseif pmethod==2
    % check for nonconstant coefficients
    [dcoeffdx,dcoeffdy] = specific_gradcoeff(x,y,size(evt,1));
    if norm(dcoeffdx)+ norm(dcoeffdx)~=0
        error('Oops.. P2 with variable coefficients is not allowed!');
    end
    [xy,evt,bound]=p2grid(xy,evt,bound,eboundt);
    [A,M,f] = femp2_diff(xy,evt,xy,evt);
end 

%% boundary conditions
[Agal,fgal] = nonzerobc(A,f,xy,bound);

%% save resulting system
fprintf('system saved in system_adiff.mat ...\n')
gohome
cd datafiles
save system_adiff pmethod Agal M  fgal  xy

%% compute solution
tic, fprintf('solving linear system ...  ')
x_gal=Agal\fgal;
fprintf('done\n')
etoc=toc; fprintf('Galerkin system solved in  %8.3e seconds\n',etoc)
save system_adiff x_gal  -append

%% plot solution
xx=xy(:,1); yy=xy(:,2);
xxmin=min(xx); xxmax=max(xx); xxd=(xxmax-xxmin)/256;
x=(xxmin:xxd:xxmax);
yymin=min(yy); yymax=max(yy); yyd=(yymax-yymin)/256;
y=(yymin:yyd:yymax);
tsolplot(pmethod,x_gal,evt,xy,x,y,12,domain);
drawnow

%% compute a posteriori error estimate
if pmethod==1
    fprintf('\n');
    pestim = default('Error estimation: linear/quadratic bubble functions 1/2? (default 1)',1);
    if pestim == 1
        % Linear midpoint hat functions
        fprintf('\n<strong>diffpost using 3 edge midpoint linear hat functions...</strong>\n');
        tic; [elerr_p,fe,ae] = diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,x_gal,subdivPar); toc;
        %fprintf('sanity check %10.4e\n',norm(elerr_p-elerr_x,inf))
        [err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
    elseif pestim == 2
        % Quadratic midpoint bubble functions
        fprintf('\n<strong>diffpost using 4 quadratic bubble functions...</strong>\n');
        tic; [elerr_p,fe,ae] = diffpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,x_gal); toc;
        %fprintf('sanity check %10.4e\n',norm(elerr_p-elerr_x,inf))
        [err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
    else
        error('Value not valid! Choose between either 1 or 2!');
    end
    % end inner if
    terrplotl(pmethod,x_gal,elerr_p,evt,xy,x,y,101);
    fprintf('<strong>Estimated energy error: %10.4e</strong>\n',norm(elerr_p,2));

elseif pmethod==2
    %% check performance with that of legacy code
    %tic,[elerr_x_p4,error_total_p4,fe_p4,ae_p4] = ...
    %    diffpost_p2_with_p4_x(xy,evt,eboundt,x_gal); toc
    tic,[elerr_p_p4,fe,ae] = ...
    diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal); toc
    %fprintf('sanity check %10.4e\n',norm(elerr_p_p4-elerr_x_p4,inf))
    terrplotl(pmethod,x_gal,elerr_p_p4,evt,xy,x,y,102);
end
% end if
