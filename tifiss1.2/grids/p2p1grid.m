function [xy,evt,xyp,bound,eboundt]=p2p1grid(xy,evt,eboundt,bound,ifplot)
%P2P1GRID   P2-P1 element grid generator
% [x,y,xy,evt,xyp,bound,eboundt]=p2p1grid(x,y,xy,mv,mbound,bound,ifplot);
%   input
%          xy         p1 nodal coordinate vector 
%          evt        p1 element mapping matrix
%          eboundt    element boundary mapping matrix
%          bound      p1 boundary node vector
%          ifplot     if the mesh needs to be plotted (Yes/No: 0/1)
%   output
%          xy       P2 nodal coordinate vector
%          evt      P2 element node mapping matrix
%          xyp      P1 pressure node coordinate vector
%          bound    P2 boundary node vector
%          eboundt  element boundary mapping matrix
%    TIFISS function: DJS; 14 January 2016.
% Copyright (c) 2009 D.J. Silvester, Qifeng Liao
fprintf('generalizing P2-P1 mesh...\n')
xyp=xy;
[xy,evt,bound]=p2grid(xy,evt,bound,eboundt,0);

% plotting of the grid 
if ifplot
    nvtx=length(xy(:,1));
   	adj=sparse(nvtx,nvtx);
    nel=length(evt(:,1));
    for i=1:nel
    adj(evt(i,1),evt(i,2)) =1;
    adj(evt(i,2),evt(i,3)) =1;
    adj(evt(i,3),evt(i,1)) =1;
    end
    figure(9)
    gplot(adj,xy,'b')
    axis('square')
    hold on
    plot(xy(:,1),xy(:,2),'ro')
    plot(xyp(:,1),xyp(:,2),'k*',xyp(:,1),xyp(:,2),'ro')
    axis('off')
    xybd=xy(bound,:);
    plot(xybd(:,1),xybd(:,2),'ko')
    hold off
    title('P2-P1 finite element subdivision');
    drawnow
end
fprintf('mesh generalized\n')
return
