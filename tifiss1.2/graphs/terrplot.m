function terrplot(qmethod,sol,eldata,evt,xy,x,y,figno,axis_value)
%TERRPLOT   plots solution and error estimate on square-shaped domain
%   terrplot(qmethod,sol,eldata,evt,xy,x,y,figno,axis_value);
%   input
%          qmethod    approximation method
%          sol        nodal solution vector
%          eldata     element error vector
%          evt        element mapping matrix
%          xy         vertex coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          fig        figure number
%          axis_value axis for plotting
%    TIFISS function: DJS; 11 March 2017.
% Copyright (c) 2011 D.J. Silvester, Qifeng Liao
if nargin<9, axis_value=[min(x),max(x),min(y),max(y)]; end

fprintf('plotting solution and estimated errors: ')
fprintf('this might be a little slow ... ')
if qmethod==2
    ev=evt(:,1:3);
else
    ev=evt;
end
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
figure(figno)
subplot(221),contour(X,Y,xysol,20),axis('square')
title(['solution with P',num2str(qmethod),' approximation'])
axis('off'),  
subplot(222),
trimesh(ev(:,1:3),xy(:,1),xy(:,2),sol),
axis(axis_value,'square')
view(330,30)
%%
xx=xy(:,1); yy=xy(:,2);
nel=length(eldata);
% loop over elements    
for ielem = 1:nel
xl = xx(ev(ielem,:));
yl = yy(ev(ielem,:)); 
xc(ielem,1) = (1/3)*sum(xl);
xc(ielem,2) = (1/3)*sum(yl);
end
%
% interpolate to a cartesian product mesh
x=0.5*(x(1:end-1)+x(2:end));
y=0.5*(y(1:end-1)+y(2:end));
[X,Y]=meshgrid(x,y);
xysol = griddata(xc(:,1),xc(:,2),eldata,X,Y);
subplot(223),contour(X,Y,xysol,50),axis('square')
title('estimated error')
subplot(224),mesh(X,Y,xysol),
axis(axis_value,'square')
view(330,30)
subplot(111)
fprintf('done\n')
return
