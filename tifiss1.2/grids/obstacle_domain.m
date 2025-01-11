function obstacle_domain
%OBSTACLE_DOMAIN obstacle domain triangular grid generator
%   input
%          square_type   1 for  [0,1]x[0,1]
%                        2 for [-1,1]x[-1,1]               
%          grid_reg      1 for uniform grid
%                        0 for option to stretch grid              
%
% Call DISTMESH function: distmesh_2d
% grid defining data is saved to the file: obstacle_grid.mat
%    TIFISS function: LR; 24 January 2019.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao

cylinder_choice=default('square/circle obstacle 1/2? (default circle)',2);

% set DISTMESH parameters
h=default('target mesh size (default 0.1)',0.1);
fprintf('generating triangulation ...\n')
if cylinder_choice==1
   [xy,evt]= distmesh2d(@square_obstacle_fd, @square_obstacle_fh, h, [0,-1;8,1],...
   [1.75,-0.25;2.25,-0.25;2.25,0.25;1.75,0.25]);
   elseif cylinder_choice==2
   [xy,evt]= distmesh2d(@circle_obstacle_fd, @circle_obstacle_fh, h, [0,-1;8,1],[]);
   else
   error('Oops. obstacle shape is not defined')
end
fprintf('done\n'), 

%% get boundary nodes and edges
[eboundt,bound]=find_boundary_elements(evt);
%% check the mesh regularity
[xy,evt,eboundt,bound] = mesh_regularity_check(xy,evt,eboundt,bound);

%% check the boundary nodes
neboundt=length(eboundt(:,1));
for i=1:neboundt
    iel=eboundt(i,1);
    nodes=evt(iel,:);
    if eboundt(i,2)==1, nodeb=evt(iel,[2,3]); xybd=xy(nodeb,:); end
    if eboundt(i,2)==2, nodeb=evt(iel,[3,1]); xybd=xy(nodeb,:); end
    if eboundt(i,2)==3, nodeb=evt(iel,[1,2]); xybd=xy(nodeb,:); end
%    h=norm(xybd(1,:)-xybd(2,:),2);
    % vertical bound
    if abs(xybd(1,2)-xybd(2,2))>abs(xybd(1,1)-xybd(2,1))
       if xybd(1,1)<0.5, xy(nodeb,1)=0;
  %     elseif xybd(1,1)<2, xy(nodeb,1)=1.75;
  %     elseif xybd(1,1)<3, xy(nodeb,1)=2.25;
       elseif xybd(1,1)>3,xy(nodeb,1)=8;
       end
    else %horizontal bound
       if xybd(1,2)<-0.5, xy(nodeb,2)=-1;
%       elseif xybd(1,2)<0, xy(nodeb,2)=-0.25;
%       elseif xybd(1,2)<0.5, xy(nodeb,2)=0.25;
       elseif xybd(1,2)>0.5,  xy(nodeb,2)=1;
       end
    end        
end
%% outflow boundary (natural)
outbc_query=1; %default('outflow boundary: natural/prescribed 1/2 (default is  natural)',1);
if outbc_query==1,
  [eboundt,bound]=natural_out_flowbc(xy,evt,eboundt,bound);
end

gohome
cd datafiles
save obstacle_grid.mat evt xy bound eboundt cylinder_choice;
clear 
return