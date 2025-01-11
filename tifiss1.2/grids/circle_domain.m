function circle_domain
%CIRCLE_DOMAIN unit circle triangular grid generator
% calls DISTMESH function: distmesh2d
% grid defining data is saved to the file: circle_grid.mat
%    TIFISS function: DJS; 12 February 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao

% set DISTMESH parameters
h=default('target mesh size (default 0.1)',0.1);
fprintf('generating triangulation ...\n')
fd=@(p) sqrt(sum(p.^2,2))-1;
[xy,evt]= distmesh2d(fd,@huniform,h,[-1,-1;1,1],[]);
fprintf('done\n')

%% get boundary nodes and edges
[eboundt,bound]=find_boundary_elements(evt);
%% check the mesh regularity
[xy,evt,eboundt,bound] = mesh_regularity_check(xy,evt,eboundt,bound);

gohome
cd datafiles
save circle_grid.mat xy evt bound eboundt 