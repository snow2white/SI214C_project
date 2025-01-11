function [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt)
%ADJUST_UNSTRUCT_MESH renumbering nodes (i.e., edges) in evt for unstructured meshes  
%
% [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt)
%
% input:
%            evt     element mapping matrix
%             xy     nodal coordinate vector 
%        eboundt     element-edge boundary matrix 
%
% output:
%            evt     element mapping matrix
%             xy     nodal coordinate vector 
%        eboundt     element-edge boundary matrix 
%
% If an unstructured mesh has been chosen, this functions changes the order of 
% the nodes of each element in evt. In pratice,it changes the local order of 
% edges. The reason for doing this is that TIFISS as well as current mesh 
% refinement routines assume that the 2nd edge of a given element is the longest
% one and this is not guaranteed by the initial unstructured mesh generated.
%
% Note that both position of nodes and connectivity *do not* change! 
% Changes are in evt rows (horizontal swaps) and boundary edges in eboundt matrix
%
% Example:
%    evt = [4 5 1; 3 4 1; 3 1 2; 2 6 3; 3 6 4; 4 6 5]
%    xy  = [0 0; 1 0; 0.7 0.4; 0.4 0.6; 0 1; 1 1];
%    eboundt = [1 1; 3 1; 4 3; 6 1];
%    [evt,~,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt)
%    plot_mesh(evt,xy,'',1,1);
%
% Function(s) called: p1grid_detail_space
%
%   TIFISS function: LR; 24 January 2019
% Copyright (c) 2018 A. Bespalov, L. Rocchi
 
  nel = size(evt,1);   % Number of elements
  
% Recover local coordinates
  xlv = zeros(nel,3);
  ylv = zeros(nel,3);
  for ivtx = 1:3
      xlv(:,ivtx) = xy( evt(:,ivtx), 1 ); 
      ylv(:,ivtx) = xy( evt(:,ivtx), 2 );
  end

% -----------------------------------------------------------------------------  
% Compute edge lengths
% -----------------------------------------------------------------------------
  els(:,1) = sqrt( (xlv(:,3) - xlv(:,2)).^2 + (ylv(:,3) - ylv(:,2)).^2 ); % lengths of all first edges
  els(:,2) = sqrt( (xlv(:,1) - xlv(:,3)).^2 + (ylv(:,1) - ylv(:,3)).^2 ); % lengths of all second edges
  els(:,3) = sqrt( (xlv(:,2) - xlv(:,1)).^2 + (ylv(:,2) - ylv(:,1)).^2 ); % lengths of all third edges
  
% -----------------------------------------------------------------------------    
% Find the longest edge per element 
% -----------------------------------------------------------------------------
% This is done by finding the position of maximum length per element:
% - lonedge contains the longest edge per element; 
  [~,lonedge] = max(els,[],2);
  
% -----------------------------------------------------------------------------    
% Changing the order of evt node per element
% -----------------------------------------------------------------------------
  evtrep = [evt, evt(:,1), evt(:,2)];
% - if the longest edge is 1 then take columns [3,4,5] of evtrep;
% - if the longest edge is 2 then take columns [1,2,3] of evtrep, i.e., no changes;
% - if the longest edge is 3 then take columns [2,3,4] of evtrep.
  evt(lonedge==1,:) = evtrep(lonedge==1,[3,4,5]);
  evt(lonedge==2,:) = evtrep(lonedge==2,[1,2,3]);
  evt(lonedge==3,:) = evtrep(lonedge==3,[2,3,4]);  
    
% % ---------------------------------------------------------------------------
%   % !DEBUG! Check again the lengths with the new evt: maximum has to be on 2nd position 
%   % Recover local coordinates
%   xlv = zeros(nel,3);  ylv = zeros(nel,3);
%   for ivtx = 1:3
%       xlv(:,ivtx) = xy( evt(:,ivtx), 1 ); 
%       ylv(:,ivtx) = xy( evt(:,ivtx), 2 );
%   end
%   els(:,1) = sqrt( (xlv(:,3) - xlv(:,2)).^2 + (ylv(:,3) - ylv(:,2)).^2 );
%   els(:,2) = sqrt( (xlv(:,1) - xlv(:,3)).^2 + (ylv(:,1) - ylv(:,3)).^2 );
%   els(:,3) = sqrt( (xlv(:,2) - xlv(:,1)).^2 + (ylv(:,2) - ylv(:,1)).^2 );
%   [~,lonedge] = max(els,[],2);
%   if find(lonedge~=2) > 0
%       error('Bad reshaping of evt matrix');
%   end
% % ---------------------------------------------------------------------------

% -----------------------------------------------------------------------------    
% Update eboundt matrix
% -----------------------------------------------------------------------------   
% Using the detail grid for space Y automatically also finds possibles
% boundary elements with two edges on the boundary
%
% Update the detail space Y 
  [evtY,~,boundY,~] = p1grid_detail_space(xy,evt);  
  [belem,bedge] = find( ismember(evtY,boundY) );
% Update eboundt matrix with right boundary edges
  neweboundt = [belem, bedge];
    
% This step is required for those problems with Neumann boundary conditions. 
% In such case, the elements where Neumann conditions are imposed are not 
% considered boundary elements, although they have one or two edges on the
% boundary (i.e. they are not included in the eboundt matrix). This is the 
% case of, e.g., example problem OBSTACLE_ADIFF.
% 
% We check that neweboundt (which includes all boundary elements regardless 
% of boundary conditions imposed) contains the same boundary elements as the 
% input eboundt. For Neumann problems, they would not match, hence we
% remove the corresponding boundary elements from neweboundt

  if ~isequal(size(eboundt,1),size(neweboundt,1))
      % it means that eboundt contains less boundary elements. Remove the 
      % elements from neweboundt which are then not supposed 
      % to be boundary elements
      neweboundt( ismember(neweboundt(:,1),eboundt(:,1)) == 0 , : ) = [];
  end
  
% Update name 
  eboundt = neweboundt;

end % end function  