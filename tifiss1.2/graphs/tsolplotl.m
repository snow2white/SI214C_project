function tsolplotx(sol,qmethod,evt,xy,x,y,fig,domain)
%TSOLPLOTL  plots nodal data on L-shaped domain for triangular meshes
%   input
%   tsolplotl(sol,qmethod,evt,xy,x,y,fig,domain)
%          qmethod    approximation method
%          sol        nodal solution vector
%          xy         nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          axis_label axis values for plotting
%          fig        figure number
%          domain     domain type
%
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester, Qifeng Liao

fprintf('plotting solution... ')
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
[II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
if nargin>6 & domain==2
  [II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
end
figure(fig)
subplot(121),contour(X,Y,xysol,20),axis('square')
title(['solution wih P',num2str(qmethod),' approximation'])
axis('off')
subplot(122),
trimesh(evt(:,1:3),xy(:,1),xy(:,2),sol)
axis('square')
view(330,30)
subplot(111)
fprintf('done\n')
return
