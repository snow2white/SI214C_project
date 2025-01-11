function terrplotobsp1(pmethod,sol,eldata,evt,xy,x,y,figno)
%TERRPLOTOBSP1 plots solution and error estimate on obstacle domain
%   terrplotobsp1(pmethod,sol,eldata,evt,xy,x,y,figno)
%   input:
%          pmethod      approximation method
%          sol          nodal solution vector
%          eldata       element error vector
%          evt          element mapping matrix
%          xy           vertex coordinate vector  
%          x            vector of x-axis interpolation points
%          y            vector of y-axis interpolation points
%          fig          figure number
%   TIFISS function: LR; 24 January 2019
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

fprintf('plotting solution and estimated errors: ');
fprintf('this might be a little slow ... ');

% adjust element-mapping matrix
if pmethod==2
    ev=evt(:,1:3);
else
    ev=evt;
end

% load type of obstacle (1 for square, 2 for circle)
load 'obstacle_grid.mat';

% interpolate to a cartesian product mesh (solution)
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);

figure(figno);
subplot(211); contour(X,Y,xysol,20); axis equal;
title(['solution with P',num2str(pmethod),' approximation']);
axis('off');

% recover local coordinates of elements
xlv = zeros(size(ev,1),3); ylv = zeros(size(ev,1),3); 
for ivtx = 1:3
    xlv(:,ivtx) = xy(ev(:,ivtx),1);
    ylv(:,ivtx) = xy(ev(:,ivtx),2);
end
% recover elemens' centroid coordinates
xyc(:,1)=(1/3)*sum(xlv,2);
xyc(:,2)=(1/3)*sum(ylv,2);

% interpolate to a cartesian product mesh (error)
x=0.5*(x(1:end-1)+x(2:end));
y=0.5*(y(1:end-1)+y(2:end));
[X,Y]=meshgrid(x,y);
xyerr=griddata(xyc(:,1),xyc(:,2),eldata,X,Y);

% delete griddata in the circle/square obstacle
if cylinder_choice==1
    xyerr((Y>-0.25) & (Y<0.25) & (X>1.75) & (X<2.25))=nan;
else
    xyerr(sqrt((Y-0.0).^2+(X-1.75).^2)<0.25) = nan;
end

subplot(212); contour(X,Y,xyerr,80); axis equal; colorbar;
title('estimated error');
hold on; obstdomain(cylinder_choice); hold off;

fprintf('done\n');
end % end function


% child function
function obstdomain(obstype)
if obstype==1 
    % plot squared obstacle
    plot([1.75,2.25],[-0.25,-0.25],'-k'); plot([2.25,2.25],[-0.25,0.25],'-k');
    plot([2.25,1.75],[0.25,0.25],'-k'); plot([1.75,1.75],[0.25,-0.25],'-k');  
else
    % plot circled obstacle
    xx = 1.75; yy = 0.0; rr = 0.25;
    th = 0:pi/50:2*pi;
    xunit = rr * cos(th) + xx;
    yunit = rr * sin(th) + yy;
    plot(xunit,yunit,'-k');
end
% plot external rectangle domain
plot([0,8],[-1,-1],'-k'); plot([8,8],[-1,1],'-k');
plot([8,0],[1,1],'-k'); plot([0,0],[1,-1],'-k');
end % end child function