function [eboundt,bound]=find_boundary_elements(evt)
%FIND_BOUNDARY_ELEMENTS constructs eboundt and bound for Dirichlet BC
%[eboundt,bound]=find_boundary_elements(evt);
%   input
%          evt            element mapping matrix
%   output
%          eboundt        element edge boundary matrix 
%          bound          boundary nodes 
%   TIFISS function: QL; 3 April 2011.
% Copyright (c) 2011 D.J. Silvester and Q. Liao


%% Exterior boundary 
eboundt=[];
bound=[];
nel=length(evt(:,1));
for i=1:nel
    %edge1
    t=find((evt(:,1)==evt(i,2)|evt(:,2)==evt(i,2)|evt(:,3)==evt(i,2))...
          &(evt(:,1)==evt(i,3)|evt(:,2)==evt(i,3)|evt(:,3)==evt(i,3)));
    if length(t)<2, eboundt=[eboundt;i,1]; bound=[bound,evt(i,2),evt(i,3)]; end
    %edge2
    t=find((evt(:,1)==evt(i,3)|evt(:,2)==evt(i,3)|evt(:,3)==evt(i,3))...
          &(evt(:,1)==evt(i,1)|evt(:,2)==evt(i,1)|evt(:,3)==evt(i,1)));
    if length(t)<2, eboundt=[eboundt;i,2]; bound=[bound,evt(i,3),evt(i,1)]; end
    %edge3
    t=find((evt(:,1)==evt(i,1)|evt(:,2)==evt(i,1)|evt(:,3)==evt(i,1))...
          &(evt(:,1)==evt(i,2)|evt(:,2)==evt(i,2)|evt(:,3)==evt(i,2)));
    if length(t)<2, eboundt=[eboundt;i,3]; bound=[bound,evt(i,1),evt(i,2)]; end
    
end

% remove multi boundary node
for i=1:length(bound)-1
    for j=i+1:length(bound)
        if bound(j)==bound(i)
            bound(j)=nan;
        end
    end
end
t=find(bound>0);
bound=bound(t)';

end

