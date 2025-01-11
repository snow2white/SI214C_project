function [elerr_primal, ederr_primal, errest_primal, ...
          elerr_dual,   ederr_dual,   errest_dual] = ...
          goafem_diffpost_p1_with_p1_2level(xy,evt,eboundt,evtY,xyY,boundY,Ybasis,ugal,zgal,subdivPar)
%GOAFEM_DIFFPOST_P1_WITH_P1_2LEVEL two-level error estimation for both primal and dual solutions
%
% [elerr_primal, ederr_primal, errest_primal, ...
%  elerr_dual,   ederr_dual,   errest_dual] = ...
%  goafem_diffpost_p1_with_p1_2level(xy,evt,eboundt,evtY,xyY,boundY,Ybasis,ugal,zgal,subdivPar)
%
% input:
%               xy    vertex coordinate vector  
%              evt    element mapping matrix
%          eboundt    element boundary mapping matrix
%             evtY    element mapping matrix for midpoints
%              xyY    vertex coordinate vector for midpoints
%           boundY    boundary midpoints vector
%           Ybasis    Y-basis element-positions matrix
%             ugal    primal P1 fe solution
%             zgal    dual P1 fe solution
%        subdivPar    red or bisec3 uniform sub-division flag
%
% output:
%     elerr_primal    vector of element indicators (primal)
%     ederr_primal    vector of edge indicators    (primal)
%    errest_primal    global error estimate        (primal)
%       elerr_dual    vector of element indicators (dual)
%       ederr_dual    vector of edge indicators    (dual)
%      errest_dual    global error estimate        (dual)
%
% Main reference for this estimator (in the stochastic framework):
% [BPRR18] Bespalov, Praetorius, Rocchi, Ruggeri, Goal-oriented error estimation 
% and adaptivity for elliptic PDEs with parametric or uncertaint inputs, 2018. 
%
% Function(s) called: goafem_diffpost_p1_with_p1_contrib
%                     goafem_specific_bc
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nvtx = size(xy,1);      % Number of vertices  (i.e., number of X-basis functions) 
  nyb  = size(Ybasis,1);  % Number of midpoints (i.e., number of Y-basis functions (boudary ones included))

% -----------------------------------------------------------------------------                 
% STEP 1: elementwise contributions from uniform red/bisec3 refinement type
% -----------------------------------------------------------------------------                 
  [ae,bde,~,fde,gde] = goafem_diffpost_p1_with_p1_contrib(xy,evt,subdivPar);  
  
% -----------------------------------------------------------------------------                 
% STEP 2: assembling ||a0^(1/2) \grad \phi_j||^2_L^2(D), for all \phi_j in Y             
% -----------------------------------------------------------------------------  
% Extract elements and associated edges for Y-basis functions
  elems = Ybasis(:,[1,2]);   
  edges = Ybasis(:,[3,4]);
   
% Indices: elements, Y-basis, Y-basis
  ind1 = {elems(:,1), edges(:,1), edges(:,1)};
  ind2 = {elems(:,2), edges(:,2), edges(:,2)};
  
% Summing contributions  
  B0phi = ae( sub2ind(size(ae) , ind1{:}) ) + ae( sub2ind(size(ae) , ind2{:}) );
% For a Y-basis on the boundary, we have to halve the contribution, since in 
% the above line it has been doubled
  B0phi(boundY) = B0phi(boundY) ./ 2;
     
% -----------------------------------------------------------------------------                       
% STEP 3: assembling F(v) and G(v) on the rhs for all basis \phi_j in V_Y
% -----------------------------------------------------------------------------
% Indices: elements, Y-basis
  ind1 = {elems(:,1), edges(:,1)};
  ind2 = {elems(:,2), edges(:,2)};
  
% Summing contributions (primal problem) 
  ff = fde( sub2ind(size(fde) , ind1{:}) ) + fde( sub2ind(size(fde) , ind2{:}) );
% For a Y-basis on the boundary, we have to halve the contribution, since in 
% the above line it has been doubled
  ff(boundY) = ff(boundY) ./ 2; 
 
% Summing contributions (dual problem) 
  gg = gde( sub2ind(size(gde) , ind1{:}) ) + gde( sub2ind(size(gde) , ind2{:}) );
% For a Y-basis on the boundary, we have to halve the contributions since in the 
% above line we doubled the contribution coming from the only boundary element
  gg(boundY) = gg(boundY) ./ 2; 
  
% ------------------------------------------------------------------------------                               
% STEP 4: assembling B(u_X,v) and B(z_X,v), for all basis \phi_j in Y
% ------------------------------------------------------------------------------ 
% The term B(u_X,v) as well as B(z_X,v) will be a vector nyb-by-1
  Km = sparse(nyb,nvtx);
  for kYb = 1:3
      indY = evtY(:,kYb);	 
      for kXb = 1:3
          indX = evt(:,kXb);
          Km = Km  + sparse(indY, indX, bde(:, kYb, kXb), nyb, nvtx);
      end
  end  
  
% Global rhs for primal (error) problem: F(v) - B(u_X,v)
  Rhs  = ff - Km * ugal;
  
% Global rhs for dual (error) problem: G(v) - B(z_X,v) 
  Goal = gg - Km * zgal;  

% -----------------------------------------------------------------------------  
% STEP 5:Impose interpolated error on all Y-boundary nodes
% -----------------------------------------------------------------------------                                
  [Rhs,Goal,~,nonzerobndY] = nonzerobc_yspace(evt,evtY,xy,xyY,eboundt,Rhs,Goal);
  
% -----------------------------------------------------------------------------                                  
% STEP 6: Compute the global 2-level estimate  
% -----------------------------------------------------------------------------                                  
% Interior midpoints (Y-nodes)
  interiorY = 1:boundY(1)-1;
  
% Contribution from internal Y-basis functions (primal and dual)
  yp_err_est_sq_primal = sum( (Rhs(interiorY).^2)  ./ B0phi(interiorY) );
  yp_err_est_sq_dual   = sum( (Goal(interiorY).^2) ./ B0phi(interiorY) );

% Contribution from the Y-boundary basis associated with nonzero bc (primal and
% dual): in case of homogeneous conditions this contributions is zero
  yp_err_est_sq_primal = yp_err_est_sq_primal + sum( (Rhs(nonzerobndY).^2)  ./ B0phi(nonzerobndY) );
  yp_err_est_sq_dual   = yp_err_est_sq_dual   + sum( (Goal(nonzerobndY).^2) ./ B0phi(nonzerobndY) );

% Total 2-level estimate (primal and dual)
  errest_primal = full( sqrt(yp_err_est_sq_primal) );
  errest_dual   = full( sqrt(yp_err_est_sq_dual)   );

% -----------------------------------------------------------------------------
% STEP 7: recovering the elementwise indicators from edge indicators
% -----------------------------------------------------------------------------  
% Vector of 2-level edge-indicators (primal and dual)
  ederr_primal   = sqrt( ( Rhs.^2 )  ./ B0phi );  
  ederr_dual     = sqrt( ( Goal.^2 ) ./ B0phi );
%
% Vector of 2-level element-indicators (primal and dual)
  [elerr_primal] = get_yp_errelem(evtY,ederr_primal);
  [elerr_dual]   = get_yp_errelem(evtY,ederr_dual);  

end % end function


% ------------------------------------------------------------------------------
% Child function
% ------------------------------------------------------------------------------
function [Rhs,Goal,zerobndY,nonzerobndY] = nonzerobc_yspace(evt,evtY,xy,xyY,eboundt,Rhs,Goal)
%Set boundary conditions on edges-midpoints (Y-space)

% Extract boundary elements with boundary edges respectively = 1, 2, and 3
  beled1 = eboundt( eboundt(:,2) == 1 , 1);
  beled2 = eboundt( eboundt(:,2) == 2 , 1);
  beled3 = eboundt( eboundt(:,2) == 3 , 1);
  
% For boundary edge = 1, the X-nodes positions = (2,3) and Y-node position = 1
  bnodesX1 = evt ( beled1, [2 3]);
  bnodeY1  = evtY( beled1, 1);
  [error1,bc_nodeY1] = goafem_interperror_ybc(xy,xyY,bnodesX1,bnodeY1);
 
% For boundary edge = 2, the X-nodes positions = (3,1) and Y-node position = 2
  bnodesX2 = evt ( beled2, [3 1]);
  bnodeY2  = evtY( beled2, 2);
  [error2,bc_nodeY2] = goafem_interperror_ybc(xy,xyY,bnodesX2,bnodeY2);

% For boundary edge = 3, the X-nodes positions = (1,2) and Y-node position = 3
  bnodesX3 = evt ( beled3, [1 2]);
  bnodeY3  = evtY( beled3, 3);
  [error3,bc_nodeY3] = goafem_interperror_ybc(xy,xyY,bnodesX3,bnodeY3);
  
% Update the vector Fv in all boundary-node positions  
  Rhs(bnodeY1) = error1;
  Rhs(bnodeY2) = error2;
  Rhs(bnodeY3) = error3;
%  
  Goal(bnodeY1) = error1;
  Goal(bnodeY2) = error2;
  Goal(bnodeY3) = error3;
    
% DEBUG: the following matrices contain on each row the Y-boundary node, the
% corresponding boundary condition, and the corresponding error. They are
% divided w.r.t. the Y-boundary nodes with positions 1,2, or 3, respectively.
%   [bnodeY1, bc_nodeY1, error1]
%   [bnodeY2, bc_nodeY2, error2]
%   [bnodeY3, bc_nodeY3, error3]
% Note that boundY = { bnodeY1, bnodeY2, bnodeY3 }

% -----------------------------------------------------------------------------
% STEP 2: separate Y-boundary nodes with nonzero boundary conditions
% -----------------------------------------------------------------------------
% Among all the Y-boundary nodes, we want to return those ones for which 
% the corresponding boundary conditions are non-homogeneous (if any).
% However, we have to ensure that where bc == 0 the returned value is really 
% zero (and not some value ~1.3e-17) because this would not be read as zero. 
% To this end, if some values is returned as ~1.3e-17, we find it by looking
% at those values smaller than, for instance, 1e-10.

% Set of Y-boundary nodes with zero boundary conditions
  zerobndY1 = bnodeY1( bc_nodeY1 <= 1e-10 );  
  zerobndY2 = bnodeY2( bc_nodeY2 <= 1e-10 );  
  zerobndY3 = bnodeY3( bc_nodeY3 <= 1e-10 );  
  zerobndY  = [zerobndY1; zerobndY2; zerobndY3];
  
% Set of Y-boundary nodes with non-zero boundary conditions  
  nonzerobndY1 = bnodeY1( ~ismember(bnodeY1,zerobndY1) );
  nonzerobndY2 = bnodeY2( ~ismember(bnodeY2,zerobndY2) );
  nonzerobndY3 = bnodeY3( ~ismember(bnodeY3,zerobndY3) );
  nonzerobndY  = [nonzerobndY1; nonzerobndY2; nonzerobndY3];
   
end % end child function


% ------------------------------------------------------------------------------
% Child function
% ------------------------------------------------------------------------------
function [error,bc_nodeY] = goafem_interperror_ybc(xy,xyY,bnodesX,bnodeY)
%Computes the error committed at boundary midpoints

% Get the coordinates of the boundary X-nodes and of the boundary Y-node
  allxbd_X = reshape( xy(bnodesX,1), size(bnodesX,1), 2);
  allybd_X = reshape( xy(bnodesX,2), size(bnodesX,1), 2);
  xybd_Y   = xyY(bnodeY, :);
  
% Compute the boundary conditions for the given nodes
  [bc_firstNodeX]  = goafem_specific_bc(allxbd_X(:,1), allybd_X(:,1));
  [bc_secondNodeX] = goafem_specific_bc(allxbd_X(:,2), allybd_X(:,2));
  [bc_nodeY]       = goafem_specific_bc(xybd_Y(:,1),   xybd_Y(:,2));
 
% Interpolated error
  error = bc_nodeY - 0.5 * (bc_firstNodeX + bc_secondNodeX);

end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [yp_err_el] = get_yp_errelem(evtY,edgevec)
%Recovers elementwise spatial indicators by square summing the 
%3 edge indicators per element
  
% Get the matrix nel-by-3 ypedgelem having per each row (element) the 
% correspoding edge-indicator (column)
  evtYvec = reshape(evtY',3*size(evtY,1),1);
  edgelem = edgevec(evtYvec);
  edgelem = reshape(edgelem,3,length(edgelem)/3)';
    
% Compute the elementwise estimates by square summing the 3 edge 
% indicators per element
  yp_err_el = sqrt( sum(edgelem.^2, 2) );
          
end % end child function 