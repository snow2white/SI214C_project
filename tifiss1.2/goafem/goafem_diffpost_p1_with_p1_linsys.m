function [elerr_primal, ederr_primal, est_primal, ...
          elerr_dual,   ederr_dual,   est_dual] = ...
          goafem_diffpost_p1_with_p1_linsys(xy,evt,eboundt,evtY,xyY,boundY,Ybasis,ugal,zgal,subdivPar)
%GOAFEM_DIFFPOST_P1_WITH_P1_LINSYS computes spatial error estimation solving the fully assembled linear system
%
% [elerr_primal, ederr_primal, est_primal, ...
%  elerr_dual,   ederr_dual,   est_dual] = ...
%  goafem_diffpost_p1_with_p1_linsys(xy,evt,eboundt,evtY,xyY,boundY,Ybasis,ugal,zgal,subdivPar)
%
% input:
%               xy     vertex coordinate vector  
%              evt     element mapping matrix
%          eboundt     element edge boundary matrix
%             evtY     element mapping matrix for midpoints
%              xyY     vertex coordinate vector for midpoints
%           boundY     boundary midpoints vector
%           Ybasis     Y-basis element-positions matrix
%             ugal     primal P1 fe solution
%             zgal     dual P1 fe solution
%        subdivPar     red or bisec3 uniform sub-division flag
%
% output:
%      elerr_primal    vector of element indicators (primal)
%      ederr_primal    vector of edge indicators    (primal)
%     errest_primal    global error estimate        (primal)
%        elerr_dual    vector of element indicators (dual)
%        ederr_dual    vector of edge indicators    (dual)
%       errest_dual    global error estimate        (dual)
%
% The function solves the discrete problems for the eY estimators:
%
%   B(eY,v) = F(v) - B(ugal,v),   for all v \in Y,        (1)
%   B(eY,v) = G(v) - B(zgal,v),   for all v \in Y.        (2)
%
% Recall that the source F(v) (resp. G(v)) of (1) (resp. of (2)) follows the 
% representation of [MS09] (and [FPZ16]); see also GOAFEM_FEMP1_ADIFF.
%
% NOTE that given the midpoint-solutions (edge solutions) to problem (1)
% (resp. (2)) per each node, the element indicators are computed by definition 
% of local B-norms restricted to single elements.
%
% References:
%
% [MS09] Mommer, Stevenson, A goal-oriented finite element method with
% convergence rates, SIAM J. Numer. Anal., 47(2)861-866, 2009;
%
% [FPZ16] Feischl, Praetorius, van der Zee, An abstract analysis of optimal 
% goal-oriented adaptivity, SIAM J. Numer. Anal., 54(3)1423-1448, 2016;
%
% Function(s) called:  goafem_diffpost_p1_with_p1_detcontrib
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  nvtx = size(xy,1);      % Number of vertices  (i.e., number of X-basis functions) 
  nyb  = size(Ybasis,1);  % Number of midpoints (i.e., number of Y-basis functions, boudary ones included)
  
% -----------------------------------------------------------------------------                 
% STEP 1: elementwise contributions from uniform red/bisec3 refinement type
% -----------------------------------------------------------------------------                 
  [ade,bde,~,fde,gde] = goafem_diffpost_p1_with_p1_contrib(xy,evt,subdivPar);

% -----------------------------------------------------------------------------
% STEP 2: assembling the LHS matrix (the same for both primal and dual problems)
% -----------------------------------------------------------------------------
  A = sparse(nyb,nyb);
  for i = 1:3
      % Indices for Y-basis 
      idxY1 = evtY(:,i);
      for j = 1:3
          % Indices for Y-basis
          idxY2 = evtY(:,j);
          % Assembling
          A = A + sparse(idxY1, idxY2, ade(:,i,j), nyb, nyb);
      end
  end
     
% -----------------------------------------------------------------------------                       
% STEP 3: assembling F(v) and G(v) on the rhs for all basis v in Y
% -----------------------------------------------------------------------------
% Extract elements and associated edges for Y-basis functions
  elems = Ybasis(:,[1,2]);   
  edges = Ybasis(:,[3,4]);
  
% Indices: elements, Y-basis.
  idx1 = {elems(:,1), edges(:,1)}; 
  idx2 = {elems(:,2), edges(:,2)};

% Summing deterministic contributions (primal problem) 
  ff = fde( sub2ind(size(fde) , idx1{:}) ) + fde( sub2ind(size(fde) , idx2{:}) );
% For a Y-basis on the boundary, we have to halve the contributions since in the 
% above line we doubled the contribution coming from the only boundary element
  ff(boundY) = ff(boundY) ./ 2;
  
% Summing deterministic contributions (dual problem) 
  gg = gde( sub2ind(size(gde) , idx1{:}) ) + gde( sub2ind(size(gde) , idx2{:}) );
% For a Y-basis on the boundary, we have to halve the contributions since in the 
% above line we doubled the contribution coming from the only boundary element
  gg(boundY) = gg(boundY) ./ 2; 
  
% ------------------------------------------------------------------------------                               
% STEP 4: assembling B(ugal,v) and B(zgal,v), for all basis v in Y
% ------------------------------------------------------------------------------ 
% The terms B(ugal,v) and B(zgal,v) will be a vector nyb-by-1
  Km = sparse(nyb,nvtx);
  for kYb = 1:3
      indY = evtY(:,kYb);	 
      for kXb = 1:3
          indX = evt(:,kXb);
          Km = Km  + sparse(indY, indX, bde(:, kYb, kXb), nyb, nvtx);
      end
  end  
  
% Global rhs for primal (error) problem: F(v) - B(ugal,v)
  Rhs  = ff - Km * ugal;
  
% Global rhs for dual (error) problem: G(v) - B(zgal,v) 
  Goal = gg - Km * zgal;
  
% -----------------------------------------------------------------------------  
% STEP 5: imposing interpolated error as Dirichlet bcs on Y-boundary nodes
% -----------------------------------------------------------------------------                                
  [Rhs,Goal] = diffpost_nonzerobc(evt,evtY,xy,xyY,boundY,eboundt,A,Rhs,Goal);

% -----------------------------------------------------------------------------
% STEP 6: solving the associated full system for interior positions
% -----------------------------------------------------------------------------
% The final ederr_primal/dual will contain the error solution at each midpoint: 
% note that they are nothing but that the edge (or midpoint) indicators.
  ederr_primal = zeros(nyb,1);
  ederr_dual   = zeros(nyb,1);

% Set of interior midpoints 
  interiorY = 1:boundY(1)-1;
   
% Solving the linear systems for interior positions (primal and dual);
% NOTE that there should be A' in the dual system, but A is symmetric;   
  ederr_primal(interiorY,:) = A(interiorY,interiorY) \ Rhs(interiorY,:);
  ederr_dual(interiorY,:)   = A(interiorY,interiorY) \ Goal(interiorY,:);

% Update the error values in boundary-midpoints position
  ederr_primal(boundY,:) = Rhs(boundY,:);
  ederr_dual(boundY,:)   = Goal(boundY,:);

% -----------------------------------------------------------------------------
% STEP 7: recovering the element indicators from edge indicators
% -----------------------------------------------------------------------------  
  [elerr_primal] = get_errelem(ade,evtY,ederr_primal);
  [elerr_dual]   = get_errelem(ade,evtY,ederr_dual);  
  
% -----------------------------------------------------------------------------  
% STEP 8: global estimates ||eY||_B0^2 = B0(eY,eY) (primal and dual)
% -----------------------------------------------------------------------------  
% Compute ||eY||_B0^2 = B0(eY,eY) for both primal and dual problems
  errPrimal_est_sq = ederr_primal' * Rhs;
  errDual_est_sq   = ederr_dual'   * Goal;
% Compute ||eY||_B0 = B0(eY,eY)^(1/2) for both primal and dual problems
  est_primal = full( sqrt( errPrimal_est_sq ) );
  est_dual   = full( sqrt( errDual_est_sq ) );

end  % end function


% -----------------------------------------------------------------------------  
% Child function
% -----------------------------------------------------------------------------     
function [Rhs,Goal] = diffpost_nonzerobc(evt,evtY,xy,xyY,boundY,eboundt,K,Rhs,Goal)
%Imposes Dirichlet boundary conditions at boundary midpoints
  
  interiorY = 1:boundY(1)-1;

% -----------------------------------------------------------------------------
% Extract boundary elements with boundary edges respectively = 1, 2, and 3
% ----------------------------------------------------------------------------- 
  beled1 = eboundt( eboundt(:,2)==1 , 1);     %nbel1 = length(beled1);
  beled2 = eboundt( eboundt(:,2)==2 , 1);     %nbel2 = length(beled2);
  beled3 = eboundt( eboundt(:,2)==3 , 1);     %nbel3 = length(beled3);
   
% -----------------------------------------------------------------------------
% Extract boundary nodes/midpoints and compute the error
% -----------------------------------------------------------------------------
  node1    = evt ( beled1, [2 3]);   % nodes of the 1st boundary edges (columns 2,3 of evt)
  midp1    = evtY( beled1, 1);       % midpoints on 1st boundary edges (column 1 of evtY)
  [error1] = interperror_bc(xy,xyY,node1,midp1);
  
  node2    = evt ( beled2, [3 1]);   % nodes of the 2nd boundary edges (columns 3,1 of evt)
  midp2    = evtY( beled2, 2);       % midpoints on 2nd boundary edges (column 2 of evtY)
  [error2] = interperror_bc(xy,xyY,node2,midp2);
  
  node3    = evt ( beled3, [1 2]);   % nodes of the 3rd boundary edges (columns 1,2 of evt)
  midp3    = evtY( beled3, 3);       % midpoints on 3rd boundary edges (column 3 of evtY)
  [error3] = interperror_bc(xy,xyY,node3,midp3);
  
% -----------------------------------------------------------------------------
% Vector of interpolated error on boundary midpoints
% -----------------------------------------------------------------------------
  miderrors = sortrows([midp1, error1; midp2, error2; midp3, error3],1);
  bcY = miderrors(:,2); 
  
% -----------------------------------------------------------------------------
% Update the Rhs in the interior-midpoint positions
% -----------------------------------------------------------------------------
% We have to do the following (primal problem):
%
%    Rhs(interiorY,:) = Rhs(interiorY,:) - K(interiorY,boundY) * bcY
%
% Note that this updated rhs for interior positions is a nIntY-by-1 vector;
% The same for the dual problem but note that the matrix K has to be K'. 
% However, in our case K=K'
  Rhs(interiorY,:)  = Rhs(interiorY,:)  - K(interiorY,boundY) * bcY;
  Goal(interiorY,:) = Goal(interiorY,:) - K(interiorY,boundY) * bcY;

% Update values in boundary-midpoints position with the imposed error
  Rhs(boundY,:)  = bcY;
  Goal(boundY,:) = bcY;
  
end % end child function
 

% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [error] = interperror_bc(xy,xyY,bnodesX,bnodeY)
%Computes the error committed at boundary midpoints

% Coordinates of the boundary nodes and of boundary edge's midpoints
  allxbd_X = reshape( xy(bnodesX,1), size(bnodesX,1), 2);
  allybd_X = reshape( xy(bnodesX,2), size(bnodesX,1), 2);
  xybd_Y   = xyY(bnodeY, :);
  
% Compute boundary conditions for the given nodes
  [bc_firstNodeX]  = goafem_specific_bc( allxbd_X(:,1), allybd_X(:,1) );
  [bc_secondNodeX] = goafem_specific_bc( allxbd_X(:,2), allybd_X(:,2) );
  [bc_nodeY]       = goafem_specific_bc( xybd_Y(:,1),   xybd_Y(:,2)   );

% Interpolated error
  error = bc_nodeY - 0.5 * (bc_firstNodeX + bc_secondNodeX);
  
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [elerr] = get_errelem(ade,evtY,edgevec)
%Recovers elementwise B_0^2 energy norms

  nel      = size(evtY,1);
  errloc   = zeros(nel,3);
  elerr_sq = zeros(nel,1);
     
% Rearranging spatial nodal values elementwise
  for imid = 1:3
      errloc(:,imid) = edgevec(evtY(:,imid));
  end 
     
% Computing elementwise B_0^2 energy norms
  for imid = 1:3          
      for jmid = 1:3
          elerr_sq(:) = elerr_sq(:) + errloc(:,imid) .* ade(:,imid,jmid) .* errloc(:,jmid);
      end
  end
  
% B_0 energy norms
  elerr = sqrt( elerr_sq );
  
end % end child function 