function tsolplotx(sol,xy,evt,fig,cmap,axis_label)
%TSOLPLOTX plots nodal data for P1 and P2 approximation
%   tsolplotx(x_gal,xy,evt,fig,cmap,axis_label);
%   input
%          x_gal       nodal solution vector
%          xy          nodal coordinate vector
%          evt         element node mapping matrix
%          fig         figure number
%          cmap        colormap (optional; default is'summer')
%          axis_label  optional axis values for plotting
%
% calls function twod_plotc (written by Jeff Borggaard)
%   TIFISS function: DJS; 2 March 2017.
% Copyright (c) 2011 D.J. Silvester, Qifeng Liao
if nargin<5, cmap='parula'; end
if nargin<6,
axis_label=[min(xy(:,1)),max(xy(:,1)),min(xy(:,2)),max(xy(:,2))];end
fprintf('plotting solution... ')
% interpolate to a cartesian product mesh
figure(fig)
subplot(121),
trimesh(evt(:,1:3),xy(:,1),xy(:,2),sol)
axis(axis_label,'square')
view(330,30)
subplot(122),
twod_plotc(fig,xy,evt,sol)
colormap(cmap)
axis(axis_label,'square'),  axis ('off')
%if all([min(x),max(x),min(y),max(y)] == [-1,1,-1,1]), squarex, end
fprintf('done\n')
return
