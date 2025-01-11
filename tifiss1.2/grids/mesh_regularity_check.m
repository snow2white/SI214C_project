function [xy,evt,eboundt,bound] = mesh_regularity_check(xy,evt,eboundt,bound)
%MESH_REGULARITY_CHECK ensures each element has at most 1 edge on boundary
%[xy,evt,eboundt,bound] = mesh_regularity_check(xy,evt,eboundt,bound);
%   input
%          xy             vertex coordinate vector  
%          evt            element mapping matrix
%          eboundt        element edge boundary matrix 
%          bound          boundary nodes 
%   ouput
%          xy             vertex coordinate vector  
%          evt            element mapping matrix
%          eboundt        element edge boundary matrix 
%          bound          boundary nodes 
%   TIFISS function: QL; 7 April 2011.
% Copyright (c) D.J. Silvester and Qifeng Liao

nv=length(xy(:,1));
nel=length(evt(:,1));
%% find the elements contain two boundary edges
n=length(eboundt(:,1));
count=0;
for i=1:n-1
    for j=i+1:n
        if eboundt(i,1)==eboundt(j,1)
         count=count+1;
         k=eboundt(i,1);
         nodes =evt(k,:);
         coords=xy(nodes,:);
         xy(nv+1,1)=1/3*sum(coords(:,1));
         xy(nv+1,2)=1/3*sum(coords(:,2));
         evt(k,3)=nv+1;
         evt(nel+1,:)=[nodes(2),nodes(3),nv+1];
         evt(nel+2,:)=[nodes(3),nodes(1),nv+1];
         nv=length(xy(:,1));
         nel=length(evt(:,1));
        end
    end
end
if count>0,[eboundt,bound]=find_boundary_elements(evt); end
