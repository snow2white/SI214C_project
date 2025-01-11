function tsolplot(qmethod,sol,evt,xy,x,y,fig,domain,axis_label)
%TSOLPLOT   plots nodal data on square-shaped domain for triangular meshes
%   tsolplot(qmethod,sol,evt,xy,x,y,fig,domain,axis_label);
%   input
%          qmethod    approximation method
%          sol        nodal solution vector
%          xy         nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          axis_label axis values for plotting
%          fig        figure number
%          domain     domain type
%          axis_label axis values for plotting
%
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester, Qifeng Liao

if nargin<9, axis_label=[min(x),max(x),min(y),max(y)];end
if nargin<8, domain=nan; end

fprintf('plotting solution... ')
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
if domain==2
  [II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
end
figure(fig)
subplot(121),contour(X,Y,xysol,20),axis('square')
title(['solution with P',num2str(qmethod),' approximation'])
axis('off')
subplot(122),
trimesh(evt(:,1:3),xy(:,1),xy(:,2),sol)
axis(axis_label,'square')
view(330,30)
fprintf('done\n')
return