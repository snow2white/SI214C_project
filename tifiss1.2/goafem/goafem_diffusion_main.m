%GOAFEM_DIFFUSION_MAIN main driver for deterministic goal-oriented FEM
%   TIFISS scriptfile: LR; 24 January 2019.
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

% -------------------------------------------------------------------------  
% Parameters
% -------------------------------------------------------------------------

% Enrichment space for hierarchical error estimator (only P1 implemented)
  pestim = 1;
  
% Red/Bisec3 error estimation 1/2? See:
% Ainsworth, Oden, A posteriori error estimation in finite element analysis, 
% Wiley, 2000 - Figure 5.2 (p. 87) for the basis functions in both cases.
  subdivPar = 2;    

% Error estimation type (for P1 approximations only)
% 1 - eY hierarchical estimator (elementwise residual problems);
% 2 - eY hierarchical estimator (assembled system for the residual problems);
% 3 - 2-level error estimator.
  estype = 3;
  
% Marking elements or edges 1/2? 
% This depends on the estimation type:
% - estype=1   -> only elements can be marked, i.e., markedgelem=1;
% - estype=2/3 -> both elements and edges can be marked, i.e., markedgelem=1/2.
  markedgelem = 2;
  if ismember(markedgelem,[1,2])
      if estype == 1 && markedgelem == 2
          error('Hierarchical elementwise eY allows marking elements only!');
      end
  else
      error('Marking type not allowed!');
  end
   
% Tolerance
  fprintf('\n');
  if sn == 1
      tolerance = default('Tolerance (default 2e-05)',2e-05); 
  elseif ismember(sn,[2,3])
      tolerance = default('Tolerance (default 6e-05)',6e-05); 
  else
      tolerance = default('Tolerance (default 6e-04)',6e-04);
  end

% Marking strategy: only Doerfler marking supported;
% Six different approaches to relate primal and dual problems
  fprintf('Available strategies for Doerfler marking:');
  fprintf('\n   1.  [FPZ16] Feischl, Praetorius, Van Der Zee');
  fprintf('\n   2.  [MS09]  Mommer, Stevenson');
  fprintf('\n   3.  [HP16]  Holst, Pollock');
  fprintf('\n   4.  [BET11] Becker, Estecahandy, Trujillo');
  fprintf('\n   5.  Only primal problem');
  fprintf('\n   6.  Only dual problem\n');
  strategy  = default('Choose marking strategy: 1/2/3/4/5/6? (default 1)',1);
  threshold = default('Threshold parameter (default 0.3)',0.3);
  
% -------------------------------------------------------------------------  
% Domain generation 
% -------------------------------------------------------------------------
  fprintf('\n<strong>Coarse grid generation for a </strong>');
  if dom_type == 1
      % Unit square domain (0,1)^2, only structured mesh
      fprintf('<strong>unit square domain</strong>\n');
      square_domain_fa(1,1);
      load square_grid;
      [evt,eboundt] = p1grid(xy,mv,bound,mbound,0);    
  elseif dom_type == 2
      % L-shaped domain
      fprintf('<strong>L-shaped domain</strong>');
      mesh_type = default('\nStructured/unstructured mesh 1/2 (default 1)',1);
      if mesh_type == 1
          ell_domain_fa;
          load ell_grid;
          [evt,eboundt] = p1grid(xy,mv,bound,mbound,0);
      elseif mesh_type == 2
          ell_domain_unstructured_fa;
          load ell_grid;
          % For unstructured meshes reorder the nodes of each element in evt
          [evt,xy,eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);
      else
          error('Invalid mesh type'); 
      end
  elseif dom_type == 3
      % Large crack domain (-1,1)^2 \ (-1,0)x{0} (crack on the left)
      fprintf('<strong>crack domain</strong>');
      crack_domain_large;
      load crack_grid;
  end
    
% In case of pointwise estimation, enter point and radius for the QoI:
  if sn == 4
      fprintf('\nInternal point (x0,y0) and radius ''r'' for the QoI:\n');
      x0 = default('   x0 (default +0.40)',0.40);
      y0 = default('   y0 (default -0.50)',-0.50);
      r  = default('   r  (default  0.15)',0.15);
      goafem_mollifier(dom_type,x0,y0,r);
  end  
  
% Plot the starting mesh
  plot_mesh(evt,xy,'Initial mesh');
  pause(1);
  
% Total nodes and interior nodes
  totalnodes = 1:size(xy,1);
  interior = totalnodes(~ismember(totalnodes,bound));
  
% Allocate memory
  maxIter       = 100;
  errors        = zeros(1,maxIter);
  errors_primal = zeros(1,maxIter);
  errors_dual   = zeros(1,maxIter);
  dof           = zeros(1,maxIter);
  nels          = zeros(1,maxIter);

% -------------------------------------------------------------------------
% Adaptive Finite Element Loop: SOLVE -> ESTIMATE -> MARK -> REFINE
% -------------------------------------------------------------------------
  refcont = 1;
  loopTime = tic;

  while true
      fprintf('\n'); fprintf(num2str(repmat('<strong>-</strong>',1,60))); fprintf('\n');
      fprintf('<strong>Iteration %d\n</strong>',refcont);
      fprintf(num2str(repmat('<strong>-</strong>',1,60))); fprintf('\n');

      fprintf('<strong>Number of elements:</strong> %d\n',size(evt,1));
      fprintf('<strong>Number of vertices:</strong> %d\n',size(xy,1));
      
      % -------------------------------------------------------------------
      % SOLVE
      % -------------------------------------------------------------------
      % Set up matrices
      fprintf('Setting up P1 diffusion matrices...');
      [A,~,F,G]   = goafem_femp1_adiff(xy,evt);
      % Boundary conditions
      [Agal,Fgal] = goafem_imposebc(A,F,xy,bound);
      [Bgal,Ggal] = goafem_imposebc(A,G,xy,bound);
      fprintf('done\n');
      %
      % Compute solution
      fprintf('<strong>Solving linear systems:</strong>');
      %
      % Primal problem solution
      fprintf('\nPrimal problem...');
      solveTime = tic;
      u_gal = Agal \ Fgal;
      fprintf('solved (%.5f sec)',toc(solveTime));    
      %
      % Dual problem solution
      fprintf('\nDual problem...');
      solveTime = tic;
      z_gal = Bgal' \ Ggal;
      fprintf('  solved (%.5f sec)\n',toc(solveTime));

      % -------------------------------------------------------------------
      % ESTIMATE
      % -------------------------------------------------------------------
      fprintf('<strong>A posteriori error estimation</strong>\n');
      % Computing edge lenghts/connections
      fprintf('Computing edge lengths/connections...');
      edgeGenTime = tic;
      [eex,tve,els] = tedgegen(xy,evt);
      fprintf('done (%.5f sec)\n',toc(edgeGenTime));
      % 
      % Computing detail space Y
      fprintf('Computing detail space Y info...');
      detailTime = tic;
      [evtY,xyY,boundY,Ybasis] = p1grid_detail_space(xy,evt);
      fprintf('done      (%.5f sec)\n',toc(detailTime));
      %
      % Compute a posteriori error estimation 
      errestTime = tic; 
      %
      if pestim == 1 
          %
          % P1-error estimation
          if estype == 1
              % Hierarchical eY estimator: elementwise residual problems
              fprintf('Hierarchical eY estimator: solving elementwise residual problems\n');
              [elerr_primal, est_primal, ...
               elerr_dual,   est_dual] = ...
               goafem_diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,u_gal,z_gal,evtY,xyY,subdivPar);            
          elseif estype == 2
              % Hierarchical eY estimator: solving the assembled linear system
              fprintf('Hierarchical eY estimator: solving assembled linear system\n');
              [elerr_primal, ederr_primal, est_primal, ...
               elerr_dual,   ederr_dual,   est_dual] = ...
               goafem_diffpost_p1_with_p1_linsys(xy,evt,eboundt,evtY,xyY,boundY,Ybasis,u_gal,z_gal,subdivPar); 
          else% estype == 3
              % 2-level error estimator
              fprintf('Two-level estimator\n');
              [elerr_primal, ederr_primal, est_primal, ...
               elerr_dual,   ederr_dual,   est_dual] = ...
               goafem_diffpost_p1_with_p1_2level(xy,evt,eboundt,evtY,xyY,boundY,Ybasis,u_gal,z_gal,subdivPar);                          
          end
     
      else% pestim == 2 
          %
          % Quadratic P2 Midpoint Bubble functions: NOT IMPLEMENTED YET
          %
      end
      fprintf('Estimation took %.5f sec\n',toc(errestTime));
      fprintf('Estimated energy error (primal):   %10.4e\n',est_primal);
      fprintf('Estimated energy error (dual):     %10.4e\n',est_dual);
      fprintf('<strong>Estimated energy error (product):  %10.4e</strong>\n',est_primal*est_dual);
      
      % Save current data: errors, dofs, and elements
      errors(refcont)        = est_primal * est_dual;
      errors_primal(refcont) = est_primal;
      errors_dual(refcont)   = est_dual;
      dof(refcont)           = size(xy,1); % =length(u_gal)=length(z_gal);
      nels(refcont)          = size(evt,1);

      % Check if the estimated error reach the tolerance
      if est_primal * est_dual <= tolerance
          fprintf('\n------------------------------------------------------------\n');
          fprintf('<strong>Tolerance reached!</strong>\n');
          fprintf('------------------------------------------------------------\n');
          break;
      end
      
      % -------------------------------------------------------------------
      % MARK
      % -------------------------------------------------------------------
      if markedgelem == 1 
          % Marking elements
          [Mset] = goafem_marking_strategy(elerr_primal,elerr_dual,strategy,threshold);
      else% markedgelem == 2
          % Marking edges
          [Mset] = goafem_marking_strategy(ederr_primal,ederr_dual,strategy,threshold);
      end
      %
      % Overall set of marked elements and edges
      [MMele,MMedge] = get_all_marked_elem(Ybasis,evtY,Mset,markedgelem);
          
      % -------------------------------------------------------------------
      % REFINE
      % -------------------------------------------------------------------
      fprintf('Mesh refinement...');
      meshRefTime = tic;
      % Mesh refinement
      [evt,xy,bound,interior,eboundt] = mesh_ref(MMele,MMedge,evt,xy,bound,evtY,xyY,boundY);
      fprintf('done (%.5f sec)\n',toc(meshRefTime));        

      % Update counter
      refcont = refcont + 1;  
  
  end % end while
  looptime = toc(loopTime);

% -------------------------------------------------------------------------
% Postprocessing
% -------------------------------------------------------------------------  
  gohome; cd datafiles;

% Resizing vectors
  errors        = errors(1:refcont);
  errors_primal = errors_primal(1:refcont);
  errors_dual   = errors_dual(1:refcont);
  dof           = dof(1:refcont);
  nels          = nels(1:refcont);
  
  fprintf('<strong>Total elapsed time:</strong>                %.3f sec\n',looptime);
  fprintf('<strong>Final estimated energy error:</strong>      %3.4e\n',errors(refcont));
  fprintf('<strong>Total number of refinements:</strong>       %d\n',refcont-1);
  fprintf('<strong>Final number of elements:</strong>          %d\n',size(evt,1));
  fprintf('<strong>Final number of total vertices:</strong>    %d\n',size(xy,1));
  fprintf('<strong>Final number of interior vertices:</strong> %d\n',length(interior));

% Plot the refinement path
  goafem_convplot;

% Plot the final mesh 
  plot_mesh(evt,xy,'Final mesh');
    
% Plot solution and errors
  fprintf('\nPlot primal and dual solutions and corresponding error estimates...');
  goafem_plot_data(dom_type,u_gal,elerr_primal,evt,xy,'Primal FE solution');
  goafem_plot_data(dom_type,z_gal,elerr_dual,evt,xy,'Dual FE solution');
  fprintf('done');
  
% Save data
  save goafem_adaptive_output.mat errors errors_primal errors_dual dof nels evt xy bound;  
  fprintf('\n-> Output data saved to: datafiles/goafem_adaptive_output.mat\n\n');  

% end scriptfile