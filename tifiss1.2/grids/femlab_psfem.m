function [xy,evt,bound,ebound]= femlab_psfem(p,e,t);
%FEMLAB_PSFEM  convert FEMLAB mesh to TIFISS format
%   [xy,evt,bound,ebound]= fem_cylinder(p,e,t);
%   input
%          p,e,t      femlab 'fem' structure    
%   output
%          xy         vertex coordinate vector  
%          evt        element mapping matrix
%          bound      boundary vertex vector
%          ebound     boundary edge matrix 
%
%    PIFISS function: DJS; 31 January 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester

% extract vertex co-ordinates from FEM mesh structure
px=p(1,:);   % vertex xy -coords
py=p(2,:); 
xy=[px',py'];
evt=(t(1:3,:))'; ee=e(1:2,:)';
nel=max(size(t));     % no. elements
nvtx=max(size(px));   % no. vertices
nedge=max(size(e));   % no. boundary edges

% get labels of boundary edges        
sss=max(e(5,:));              % no. of boundary segments
b_labels=sort(e(1:2,:)',2);
xbound=sort([b_labels(:,1);b_labels(:,2)]);
bound=unique(xbound);

%% sort out element edge 
ebound=zeros(nedge,2);
% first edge
edges_1=t(2:3,:)';  
adj_1=sparse(edges_1(:,1), edges_1(:,2), [1:nel]', nvtx, nvtx);
adj_1=adj_1+adj_1';
ebd=zeros(nedge,1); for ctr=1:nedge, ebd(ctr,1)=adj_1(ee(ctr,1),ee(ctr,2));end
ebound(:,1)=ebound(:,1)+ebd;
ek=find(ebd>0); ebound(ek,2)=1;
% second edge
edges_2=t([3,1],:)';
adj_2=sparse(edges_2(:,1), edges_2(:,2), [1:nel]', nvtx, nvtx);
adj_2=adj_2+adj_2';
ebd=zeros(nedge,1); for ctr=1:nedge, ebd(ctr,1)=adj_2(ee(ctr,1),ee(ctr,2));end
ebound(:,1)=ebound(:,1)+ebd;
ek=find(ebd>0); ebound(ek,2)=2;
% third edge
edges_3=t(1:2,:)';  
adj_3=sparse(edges_3(:,1), edges_3(:,2), [1:nel]', nvtx, nvtx);
adj_3=adj_3+adj_3';
ebd=zeros(nedge,1); for ctr=1:nedge, ebd(ctr,1)=adj_3(ee(ctr,1),ee(ctr,2));end
ebound(:,1)=ebound(:,1)+ebd;
ek=find(ebd>0); ebound(ek,2)=3;
return