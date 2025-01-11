function imeshplot(xy,evt,bound,eboundt,in_bound,in_eboundt)
%IMESHPLOT   triangular mesh verification
%   imeshplot(xy,evt,bound,eboundt);
%    PIFISS function: QL; 6 March 2011. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester

fprintf('\nMesh logistics ..\n')
nvtx=length(xy(:,1)); 
fprintf('  %g nodes \n',nvtx)
nelement=length(evt(:,1));
fprintf('  %g elements \n',nelement)
nboundvtx=length(bound);
fprintf('  %g nodes on Dirichlet boundary \n',nboundvtx)
nboundedge=length(eboundt(:,1));
fprintf('  %g element edges on Dirichlet boundary \n',nboundedge)
if nargin > 4
ninboundvtx=length(in_bound);
fprintf('  %g nodes on interior interface \n\n',ninboundvtx)
end
%
%
adj=sparse(nvtx,nvtx); adx=sparse(nvtx,nvtx);
for i=1:nelement
   adj(evt(i,1),evt(i,2)) =1;
   adj(evt(i,2),evt(i,3)) =1;  
   adj(evt(i,3),evt(i,1)) =1;
end
%
%% define element edges
adjb=sparse(nvtx,nvtx);
k1=find(eboundt(:,2)==1)';
for k=eboundt(k1)
   adjb(evt(k,2),evt(k,3))=1;
end
k2=find(eboundt(:,2)==2)';
for k=eboundt(k2)
   adjb(evt(k,3),evt(k,1))=1;
end
k3=find(eboundt(:,2)==3)';
for k=eboundt(k3)
   adjb(evt(k,1),evt(k,2))=1;
end

%
figure(1)
gplot(adj,xy,'b')
hold on
stnode=int2str([1:nvtx]');
text(xy(:,1),xy(:,2),stnode)
axis('equal'),axis('off')
title('Indices of nodes of the element mesh')
axis('square')
hold off
figure(2)
gplot(adj,xy,'b')
hold on
gplot(adjb,xy,'r')
xybd=xy(bound,:);
stbd=int2str([1:nboundvtx]');
text(xybd(:,1),xybd(:,2),stbd,'color','black')
title('Indices of nodes on the Dirichlet boundary')
axis('equal'),axis('off')
hold off
if nargin>4
%% define initerior interface 
xybdi=xy(in_bound,:);
adjbi=sparse(nvtx,nvtx);
k1=find(in_eboundt(:,2)==1)';
for k=in_eboundt(k1)
   adjbi(evt(k,2),evt(k,3))=1;
end
k2=find(in_eboundt(:,2)==2)';
for k=in_eboundt(k2)
   adjbi(evt(k,3),evt(k,1))=1;
end
k3=find(in_eboundt(:,2)==3)';
for k=in_eboundt(k3)
   adjbi(evt(k,1),evt(k,2))=1;
end
figure(3)
gplot(adj,xy,'b')
hold on
gplot(adjb,xy,'r')
%xybd=xy(bound,:);
%stbd=int2str([1:nboundvtx]');
gplot(adjbi,xy,'green')
xybdi=xy(in_bound,:);
stbdi=int2str([1:ninboundvtx]');
text(xybdi(:,1),xybdi(:,2),stbdi,'color','black')
title('Indices of nodes on the Interiori Interfaces')
axis('equal'),axis('off')
hold off
end
return
