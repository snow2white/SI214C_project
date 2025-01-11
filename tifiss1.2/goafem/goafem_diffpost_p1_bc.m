function [ae,fe] = goafem_diffpost_p1_bc(ae,fe,xy,evt,eboundt,evtY,xyY)
%GOAFEM_DIFFPOST_P1_BC imposes Dirichlet bcs at boundary elements for the residual problem
%
% [ae,fe] = goafem_diffpost_p1_bc(ae,fe,xy,evt,eboundt,evtY,xyY)
%
% input:
%          ae    elementwise lhs matrices (without bc imposed)
%          fe    elementwise rhs vectors  (without bc imposed)
%          xy    vertex coordinate vector  
%         evt    element mapping matrix
%     eboundt    element edge boundary matrix 
%        evtY    element mapping matrix for midpoints
%         xyY    vertex coordinate vector for midpoints
%
% output:
%          ae    elementwise lhs matrices (with bc imposed)
%          fe    elementwise rhs vectors  (with bc imposed) 
%
% Function(s) called: specific_bc
% 
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri
    
% -----------------------------------------------------------------------------
% Extract boundary elements with boundary edges respectively = 1, 2, and 3
% -----------------------------------------------------------------------------
  beled1 = eboundt( eboundt(:,2)==1 , 1);
  beled2 = eboundt( eboundt(:,2)==2 , 1);
  beled3 = eboundt( eboundt(:,2)==3 , 1);
   
% -----------------------------------------------------------------------------
% Extract boundary nodes/midpoints and compute the error
% -----------------------------------------------------------------------------
  node1    = evt ( beled1, [2 3]);    % nodes of the 1st boundary edges (columns 2,3 of evt)
  midp1    = evtY( beled1, 1);        % midpoints on 1st boundary edges (column 1 of evtY)
  [error1] = interperror_bc(xy,xyY,node1,midp1);
  
  node2    = evt ( beled2, [3 1]);    % nodes of the 2nd boundary edges (columns 3,1 of evt)
  midp2    = evtY( beled2, 2);        % midpoints on 2nd boundary edges (column 2 of evtY)
  [error2] = interperror_bc(xy,xyY,node2,midp2);
  
  node3    = evt ( beled3, [1 2]);    % nodes of the 3rd boundary edges (columns 1,2 of evt)
  midp3    = evtY( beled3, 3);        % midpoints on 3rd boundary edges (column 3 of evtY)
  [error3] = interperror_bc(xy,xyY,node3,midp3);

% -----------------------------------------------------------------------------  
% Update the RHS in internal boundary-element midpoint positions
% -----------------------------------------------------------------------------
  fe(beled1,[2,3]) = fe(beled1,[2,3]) - ae(beled1,[2,3],1) .* repmat(error1,1,2);
  fe(beled2,[1,3]) = fe(beled2,[1,3]) - ae(beled2,[1,3],2) .* repmat(error2,1,2);
  fe(beled3,[1,2]) = fe(beled3,[1,2]) - ae(beled3,[1,2],3) .* repmat(error3,1,2);
  
% -----------------------------------------------------------------------------
% Impose Dirichlet condition on LHS matrix
% -----------------------------------------------------------------------------
% Boundary elements with boundary edge = 1 
  ae(beled1,:,1) = 0.0;
  ae(beled1,1,:) = 0.0;
  ae(beled1,1,1) = 1.0;

% Boundary elements with boundary edge = 2
  ae(beled2,:,2) = 0.0;
  ae(beled2,2,:) = 0.0;
  ae(beled2,2,2) = 1.0;
  
% Boundary elements with boundary edge = 3
  ae(beled3,:,3) = 0.0;
  ae(beled3,3,:) = 0.0;
  ae(beled3,3,3) = 1.0;
     
% -----------------------------------------------------------------------------  
% Finally update the RHS in all boundary-midpoint positions
% -----------------------------------------------------------------------------
  fe(beled1,1) = error1;
  fe(beled2,2) = error2;
  fe(beled3,3) = error3;
  
end  % end function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [error] = interperror_bc(xy,xyY,bnodesX,bnodeY)
% Impose the error on boundary midpoints

% Coordinates of the boundary nodes and of boundary edge's midpoints
  allxbd_X = reshape( xy(bnodesX,1), size(bnodesX,1), 2);
  allybd_X = reshape( xy(bnodesX,2), size(bnodesX,1), 2);
  xybd_Y   = xyY(bnodeY, :);
  
% Compute boundary conditions for the given nodes
  [bc_firstNodeX]  = goafem_specific_bc(allxbd_X(:,1), allybd_X(:,1));
  [bc_secondNodeX] = goafem_specific_bc(allxbd_X(:,2), allybd_X(:,2));
  [bc_nodeY]       = goafem_specific_bc(xybd_Y(:,1),   xybd_Y(:,2));
 
% Interpolated error
  error = bc_nodeY - 0.5 * (bc_firstNodeX + bc_secondNodeX);
    
end % end child function