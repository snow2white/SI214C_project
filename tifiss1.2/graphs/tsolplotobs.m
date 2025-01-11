function tsolplotobs(qmethod,sol,evt,xy,x,y,fig,axis_label)
%TSOLPLOTOBS   plots nodal data on obstacle domain for triangular meshes
%   tsolplotobs(qmethod,sol,evt,xy,x,y,fig,axis_label);
%   input
%          qmethod    approximation method
%          sol        nodal solution vector
%          xy         nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          axis_label axis values for plotting
%          fig        figure number
%          axis_label axis values for plotting
%
%    TIFISS function: DJS; 16 January 2016.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

if nargin<8, axis_label=[min(x),max(x),min(y),max(y)];end
fprintf('plotting solution... ')
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
figure(fig)
subplot(211),contour(X,Y,xysol,20),axis('equal')
title(['solution with P',num2str(qmethod),' approximation'])
axis('off')
subplot(212),
trimesh(evt(:,1:3),xy(:,1),xy(:,2),sol)
axis(axis_label,'equal')
view(330,30)
fprintf('done\n')
return