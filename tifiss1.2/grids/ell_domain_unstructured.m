function ell_domain_unstructured
%ELL_DOMAIN_UNSTRUCTURED  L-shape unstructured grid generator
%   ell_domain_unstructured
%   input from keyboard by the default function
%    h                        the mesh size
%    mesh_density_ratio       the relative mesh density close the corner,
%                             smaller mesh_density_ratio leads to finer
%                             grids close corner
%
% calls DISTMESH function: distmesh_2d
% grid defining data is saved to the file: ell_grid.mat
%    TIFISS function: DJS; 12 February 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao

global mesh_density_ratio   % called in the function ell_fh
% set DISTMESH parameters
h=default('target mesh size (default 0.02)',0.02);
mesh_density_ratio = default('mesh density near corner (0.2)',0.2);
%it_number = 300;  % upper limit on number of mesh iterations
fprintf('generating triangulation ...\n')
[xy,evt]= distmesh2d(@ell_fd, @ell_fh, h, [-1,-1;1,1],[]);
fprintf('done\n')

%% get boundary nodes and edges
[eboundt,bound]=find_boundary_elements(evt);
%% check the mesh regularity
[xy,evt,eboundt,bound] = mesh_regularity_check(xy,evt,eboundt,bound);

%% reset boundary nodes
neboundt=length(eboundt(:,1));
for i=1:neboundt
    iel=eboundt(i,1);
    nodes=evt(iel,:);
    if eboundt(i,2)==1, nodeb=evt(iel,[2,3]); xybd=xy(nodeb,:); end
    if eboundt(i,2)==2, nodeb=evt(iel,[3,1]); xybd=xy(nodeb,:); end
    if eboundt(i,2)==3, nodeb=evt(iel,[1,2]); xybd=xy(nodeb,:); end
    h=norm(xybd(1,:)-xybd(2,:),2);
    % vertical bound or horizontal bound
    if abs(xybd(1,2)-xybd(2,2))>abs(xybd(1,1)-xybd(2,1))
       if xybd(1,1)<-0.5, xy(nodeb,1)=-1;
       elseif xybd(1,1)<0.5, xy(nodeb,1)=0;
       else xy(nodeb,1)=1;
       end
    else
       if xybd(1,2)<-0.5, xy(nodeb,2)=-1;
       elseif xybd(1,2)<0.5, xy(nodeb,2)=0;
       else  xy(nodeb,2)=1;
       end
    end        
end

x = xy(:,1); y = xy(:,2);

gohome
cd datafiles
save ell_grid.mat xy evt bound eboundt x y