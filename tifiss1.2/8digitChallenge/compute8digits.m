%COMPUTE8DIGITS adaptively solve 8-digit challenge problem
%   TIFISS scriptfile: DJS; 06 January 2019, 16 August 2019
% Copyright (c) 2018 D.J. Silvester, A. Bespalov

  gohome; cd datafiles;
% set output flag
if ~exist('outflag','var'), outflag=1; end
% set tolerance
if ~exist('tolerance','var'), tolerance = 4e-5;  end

% -------------------------------------------------------------------------  
% Parameters
% -------------------------------------------------------------------------

% fix P2 Galerkin approximation and the tolerance
pmethod = 2;
      [dcoeffdx,dcoeffdy] = specific_gradcoeff(0:0.1:1, 0:0.1:1, 50);
      if norm(dcoeffdx) + norm(dcoeffdx) ~= 0
          error('Oops .. P2 approximations with nonconstant coefficients is not allowed!');
      end

% Red/Bisec3 for spatial error estimation 1/2? See:
% Ainsworth, Oden, A posteriori error estimation in finite element analysis, 
% Wiley, 2000 - Figure 5.2 (p. 87) for the basis functions in both cases.
subdivPar = 2;

% Marking elements
fprintf('Using hierarchical estimator (elementwise residual problems)\n');
markedgelem = 1;



% Marking threshold parameters for both elements/edges and indices:
% 1 - maximum strategy:
%     large threshold -> small set of marked elements/edges
%     small threshold -> large set of marked elements/edges
% 2 - Doerfler (equilibration) strategy:
%     large threshold -> large set of marked elements/edges
%     small threshold -> small set of marked elements/edges
markstrat = 2;
threshold = 0.5;
   
% -------------------------------------------------------------------------  
% Domain generation 
% -------------------------------------------------------------------------
  fprintf('\nCoarse grid generation for a ');
  if dom_type == 1
      % 
      % Unit or reference square domain; structured mesh
      % 
      if ismember(sn,[1,3])
          % Unit square, uniform grids
          fprintf('unit square domain\n');
          square_domain_fa(1,1);
      elseif sn == 2
          % Reference grid, allowing stretched grids
          fprintf('reference square domain\n');
          square_domain_fa(1,0);
      end    
      load square_grid;
      [evt,eboundt] = p1grid(xy,mv,bound,mbound,0);
  
  elseif dom_type == 2
      % 
      % L-shaped domain (-1,1)^2 \ (-1,0]^2
      % 
      fprintf('L-shaped domain\n');
      mesh_type = 1; %structured mesh
                                 
      if mesh_type == 1
          ell_domain_fa;
          load ell_grid;
          [evt,eboundt] = p1grid(xy,mv,bound,mbound,0);
      elseif mesh_type == 2
          ell_domain_unstructured_fa;
          load ell_grid;
          % For unstructured meshes reorder the nodes of each element in evt 
          [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy);
      else
          error('Invalid mesh type'); 
      end
     
  elseif dom_type == 3
      % 
      % Crack domain (-1,1)^2 \ (-1,0)x{0} (crack on the left)
      % 
      fprintf('crack domain');
      crack_domain_large;
      load crack_grid;
  end

% Plot the starting mesh
  plot_mesh(evt,xy,'Initial mesh'); 
  pause(1);
 
% Total nodes and interior nodes
  totalnodes = 1:size(xy,1);
  interior   = totalnodes(~ismember(totalnodes,bound));
                     
% Allocate memory
  maxiter     = 30;
  comperror   = zeros(1,maxiter);       % estimated error
  energy      = zeros(1,maxiter);       % energy of the Galerkin solutions
  totdofs     = zeros(1,maxiter);       % total dofs
  nnels       = zeros(1,maxiter);       % total number of elements
  intdofs     = zeros(1,maxiter);       % internal dofs
  nnzerobdofs = zeros(1,maxiter);       % boundary dofs (with nonzero bc) 

% -----------------------------------------------------------------------------
% Adaptive Finite Element Loop: SOLVE -> ESTIMATE -> MARK -> REFINE
% -----------------------------------------------------------------------------
  refcont   = 1;
  fprintf('Stopping tolerance is %7.3e\n',tolerance); pause(1)
  totalTime = tic;
  while true
      if outflag == 1
      fprintf('\n'); fprintf(num2str(repmat('-',1,60))); fprintf('\n');
      fprintf('Iteration %d\n',refcont);
      fprintf(num2str(repmat('-',1,60))); fprintf('\n');

      fprintf('Data mesh:');
      fprintf('\nNumber of elements: %d',size(evt,1));
      fprintf('\nNumber of vertices: %d',size(xy,1));
    else
      fprintf('\n[%g]. ',refcont); end
    
      % -------------------------------------------------------------------
      % SOLVE
      % -------------------------------------------------------------------
      if outflag == 1
      fprintf('\nComputing P%d Galerkin solution\n',pmethod); end
      if pmethod == 1
          [A,M,f] = femp1_adiff(xy,evt);
      elseif pmethod == 2
          % Save oldest xy and bound data
          xyp1 = xy;
          boundp1 = bound;             
          [xy,evt,bound] = p2_grid_generator(xy,evt,bound);
          [A,M,f] = femp2_diff(xy,evt,xy,evt);
      end
      %
      % Boundary conditions
      [Agal,fgal] = nonzerobc(A,f,xy,bound);
      %
      % Compute solution
      if outflag == 1
      solveTime = tic;
      fprintf('Solving linear system...'); end
      x_gal = Agal\fgal;
      if outflag == 1
      fprintf('done (%.5f sec)\n',toc(solveTime)); end
      

%% Finding the value of x_gal at a given mesh node
%%      x0sol = zeros(1,maxiter);
%      x0ind = find(xy(:,1)==.009765625000000 & xy(:,2)==.009765625000000); %.011718750000000
%      if isempty(x0ind)==1
%         x0sol(refcont) = 0;
%      else
%         x0sol(refcont) = x_gal(x0ind);
%      end
      

% % Finding the value of u_h at a given point in the domain
       xp = 1.0e-02; yp = 1.0e-2;
       % find the element which the input point belongs to
       [el] = find_elem_given_point(xp, yp, xy(1:size(xyp1,1),:) , evt(:,1:3));
       
       % get P2 coordinates of the elements      
       coordelP2 = xy( evt(el,:), : );
       locx = coordelP2(:,1);
       locy = coordelP2(:,2);
 
       % Get P2 basis functions on the element
       BB = [ones(6,1), locx, locy, locx.^2, locx.*locy, locy.^2];
       psi = zeros(6,1);
       for k = 1:6
           b = zeros(6,1);
           b(k) = 1.0;
           
           % coefficients of the psi(k) P2 basis function on the element el
           a = BB \ b;
           
           % psi(k) evaluated at given input point
           pvalues = [1; xp; yp; xp^2; xp*yp; yp^2];
           psi(k) = a' * pvalues;
       end
       
       % then P2 x_gal solution at the given point
       P2valpoint(refcont) = x_gal( evt(el,:) )' * psi;
%

      % -------------------------------------------------------------------
      % ESTIMATE
      % -------------------------------------------------------------------
      if outflag == 1
      fprintf('A posteriori error estimation\n');
      %
      % Computing edge lenghts/connections
      fprintf('Computing edge lengths/connections...');
      edgeGenTime = tic; end
      [eex,tve,els] = tedgegen(xy,evt);
      if outflag == 1
      fprintf('done (%.5f sec)\n',toc(edgeGenTime)); end
      % 
      % Computing detail space Y
       if outflag == 1
      fprintf('Computing detail space Y info...');
      detailTime = tic; end
      [evtY,xyY,boundY,Ybasis] = p1grid_detail_space(xy,evt);
       if outflag == 1
      fprintf('done      (%.5f sec)\n',toc(detailTime));
      %
      % Compute a posteriori error estimation 
      errorTime = tic; end
      %
      if pmethod == 1 
          % P1-error estimation
          %
          if pestim == 1
              % Using linear midpoint Hat functions
              fprintf('Error estimation using 3 edge midpoint linear functions\n');
              %
              if estimtype == 1
                  % Hierarchical eY estimator: elementwise residual problem
                  fprintf('Hierarchical eY estimator: solving elementwise residual problems\n');
                  [elerr,fe,ae] = diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,x_gal,subdivPar); 
                  [err_p,elerr] = diffpost_p1_bc(ae,fe,elerr,xy,evt,eboundt);
                  %
                  % Global error estimate
                  errest = norm(elerr,2);                
              elseif estimtype == 2
                  % Hierarchical eY estimator: solving the assembled linear system
                  fprintf('Hierarchical eY estimator: solving assembled linear system\n');
                  [elerr,ederr,errest] = diffpost_p1_with_p1_linsys(evt,xy,eboundt,x_gal,evtY,xyY,boundY,Ybasis,subdivPar); 
              else% estimtype == 3
                  % 2-level error estimator
                  fprintf('Two-level estimator\n');
                  [elerr,ederr,errest] = diffpost_p1_with_p1_2level(evt,xy,eboundt,x_gal,evtY,xyY,boundY,Ybasis,subdivPar);                  
              end
    
          else%if pestim == 2  
              % Using quadratic Midpoint Bubble functions
              fprintf('Error estimation using 4 quadratic bubble functions\n');
              [elerr,fe,ae] = diffpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,x_gal);    
              [err_p,elerr] = diffpost_p1_bc(ae,fe,elerr,xy,evt,eboundt);
              %
              % Global error estimate
              errest = norm(elerr,2);
          end
    
      else% pmethod == 2
          % P2-error estimation
           if outflag == 1
          fprintf('Error estimation using quartic bubble functions\n'); end
          %
          % check performance with that of legacy code
          %tic; [elerr_x_p4,error_total_p4,fe_p4,ae_p4] = diffpost_p2_with_p4_x(xy,evt,eboundt,x_gal); toc
          [elerr,fe,ae] = diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal);
          %
          % Global error estimate
          errest = norm(elerr,2);
      end % end if
      %
       if outflag == 1
      fprintf('Estimated energy error: %10.4e\n',errest);
      fprintf('Estimation took %.5f sec\n',toc(errorTime));  end

      
      % -------------------------------------------------------------------
      % Save current iteraton data
      % -------------------------------------------------------------------
      comperror(refcont) = errest; 
      energy(refcont)    = sqrt( x_gal' * A * x_gal );   %sqrt( x_gal' * Agal * x_gal );
      nnels(refcont)     = size(evt,1);
      totdofs(refcont)   = length(x_gal);
      intdofs(refcont)   = length(x_gal) - length(bound);
      %
      % For non-homogeneous bcs (problem 5), we find the boundary nodes with 
      % nonzero bcs
      [bc] = specific_bc( xy(bound,1) , xy(bound,2) );
      zerobound    = bound(bc <= 1e-10);
      nonzerobound = bound(~ismember(bound,zerobound));
      % NOTE that bound = {zerobound, nonzerobound} and in case of homogeneous 
      % bcs, the vector nonzerobound is empty
      nnzerobdofs(refcont) = length(nonzerobound);
      

      % -------------------------------------------------------------------
      % Check if the estimated error reached the tolerance
      % -------------------------------------------------------------------
      if errest <= tolerance
          fprintf('\n------------------------------------------------------------\n');
          fprintf('Tolerance reached!\n');
          fprintf('------------------------------------------------------------\n');
          break;
      end
      
      % Fixing data if P2-approximations have been used
      if pmethod == 2
          xy    = xyp1; 
          evt   = evt(:,1:3);     
          bound = boundp1;
      end
    
      % -------------------------------------------------------------------
      % MARK
      % -------------------------------------------------------------------
      if markedgelem == 1 
          % Marking elements
          [Mset] = marking_strategy_fa(elerr,markstrat,threshold);     
      else%markedgelem == 2         
          % Marking edges
          [Mset] = marking_strategy_fa(ederr,markstrat,threshold);
      end
      %
      % Overall set of marked elements (and edges)
      [MMele,MMedge] = get_all_marked_elem(Ybasis,evtY,Mset,markedgelem);

      % -------------------------------------------------------------------
      % REFINE
      % -------------------------------------------------------------------
      if outflag == 1
      fprintf('Mesh refinement...');
      meshRefTime = tic; end
      % Mesh refinement
      [evt,xy,bound,interior,eboundt] = mesh_ref(MMele,MMedge,evt,xy,bound,evtY,xyY,boundY);
      if outflag == 1
      fprintf('done (%.5f sec)\n',toc(meshRefTime)); end

      % Update counter
      refcont = refcont + 1;  
      
  end
% end while
  looptime = toc(totalTime);

% -------------------------------------------------------------------------
% Postprocessing
% -------------------------------------------------------------------------  
  gohome; cd datafiles;
   
% Resize vectors
  comperror(refcont+1:end)   = [];      
  energy(refcont+1:end)      = [];
  totdofs(refcont+1:end)     = [];      
  nnels(refcont+1:end)       = []; 
  intdofs(refcont+1:end)     = [];
  nnzerobdofs(refcont+1:end) = []; 

% Display final data
  fprintf('Total elapsed time:                %.5f sec\n',looptime);
  fprintf('Final estimated energy error:      %3.4e\n',errest(end));
  fprintf('Total number of iterations:        %d\n',refcont);
  fprintf('Total number of refinements:       %d\n',refcont-1);
  fprintf('Final number of elements:          %d\n',nnels(end));
  fprintf('Final number of total vertices:    %d\n',totdofs(end));
  fprintf('Final number of interior vertices: %d\n',intdofs(end));
  
% Plot the final mesh
  if pmethod == 2
      % Data for "P1" (xyp1 and evt(:,1:3)) have to be used
      plot_mesh(evt(:,1:3),xyp1,'Final mesh');
  else 
      plot_mesh(evt,xy,'Final mesh');
  end
                                   
% Plot the refinement path
  convplot;
   
% Plot solution and error estimate
  fprintf('\nPlotting solution and the error estimator...');
  plot_data(pmethod,dom_type,x_gal,elerr,evt,xy);
  fprintf('done');
  
% Save data
  save adaptive_output.mat comperror energy totdofs intdofs nnzerobdofs nnels evt xy bound;  
  fprintf('\n-> Output data saved to: datafiles/adaptive_output.mat\n\n');
                                    
figure(2)
% end scriptfile
