function  [eboundt,bound]=natural_out_flowbc(xy,evt,eboundt,bound)
%NATURAL_OUT_FLOWBC treats the right boundary nodes to be interior nodes
% [eboundt,bound]=natural_out_flowbc(xy,evt,eboundt,bound)
%  input
%          xy         vertex coordinate vector  
%          evt        element vertex mapping matrix
%          evt         element mapping matrix
%          eboundt     element boundary mapping matrix
%  output
%          evt         element mapping matrix
%          eboundt     element boundary mapping matrix
%    TIFISS function: QL; 10 April 2011. 
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao

n=length(eboundt(:,1));
xmax=max(xy(:,1));
xmin=min(xy(:,1));
rbound=bound;
rnode=[];
tol=1;
for i=1:n
    t=evt(eboundt(i,1),:);
    if eboundt(i,2)==1
       enode=[t(2),t(3)];
    elseif eboundt(i,2)==2
       enode=[t(3),t(1)]; 
    else
       enode=[t(1),t(2)];
    end
    xmt=1/3*sum(xy(t,1));
    diffx=abs(xy(enode(1),1)-xy(enode(2),1));
    diffy=abs(xy(enode(1),2)-xy(enode(2),2));
    if diffx<diffy & xmt>xmax-tol
       eboundt(i,1)=nan;
       tb=find(bound==enode(1) | bound==enode(2));
       bound(tb)=nan;
       rnode=[rnode;enode(1);enode(2)];
    end
end
t=find(eboundt(:,1)>0); eboundt=eboundt(t,:);
%% fix the nodes on corners
% remove multi nodes in rnode
n=length(rnode);
for i=1:n-1
    for j=i+1:n
        if rnode(i)==rnode(j)
           rnode(j)=nan;
        end
    end
end
t=find(rnode>0); rnode=rnode(t);
rymax=max(xy(rnode,2));
rymin=min(xy(rnode,2));
t=find(xy(rnode,2)==rymax|xy(rnode,2)==rymin); corner_node=rnode(t);
for i=1:length(corner_node)
    t=find(rbound==corner_node(i)); bound(t)=rbound(t);
end
t=find(bound(:)>0);  bound=bound(t);
return