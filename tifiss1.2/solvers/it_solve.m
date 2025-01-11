%IT_SOLVE iterative solution of predefined steady problem
%    TIFISS scriptfile: DJS 13 January 2016.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

% Declare global variables for scalar problems
global amg_grid amg_smoother

if exist('pde','var')==0,
   error('You need to set up a specific discrete problem first!'), 
end
%
%
if pde==1,
   %%% DIFFUSION PROBLEM
   fprintf('discrete diffusion system ...\n')

   % select Krylov subspace method
   itmeth = default('PCG/MINRES? 1/2 (default PCG)',1);

   % set parameters
   tol = default('tolerance? (default 1e-9)',1e-9);
   maxit = default('maximum number of iterations? (default 200)',200);
   
   % select preconditioner and construct it
   fprintf('preconditioner:\n');
   fprintf('   0  none\n');
   fprintf('   1  diagonal\n');
   fprintf('   2  incomplete cholesky\n');
   fprintf('   3  algebraic multigrid\n');
   precon = default('default is AMG ',3);
   if precon==0,     % none
      M1=[]; M2=[];
   elseif precon==1, % diagonal
      D=diag(diag(Agal)); M1=sqrt(D); M2=M1; 
   elseif precon==2, % incomplete Cholesky
     M1 = ichol(Agal); M2=M1';
%    M2 = cholinc(Agal,'0'); M1=M2'; fprintf('--- cholinc.m\n')%
   elseif precon==3, % AMG
 % uses global variables amg_grid amg_smoother
      amg_grid = amg_grids_setup(Agal);
      fprintf('\nsetup done.\n')
      plot_mg = default('plot AMG grid sequence? yes/no 1/2 (default no)',2);
      if plot_mg==1, amg_coarsen_plot(amg_grid, xy); end
         smoothopt = default('PDJ/PGS smoother? 1/2 (point damped Jacobi)',1);
         if smoothopt==1
            fprintf('point damped Jacobi smoothing ..\n')
            smoother_params = amg_smoother_params(amg_grid, 'PDJ', 2);
         else
            fprintf('point Gauss-Seidel smoothing ..\n')
            smoother_params = amg_smoother_params(amg_grid, 'PGS', 2);
         end
      amg_smoother = amg_smoother_setup(amg_grid, smoother_params); 
   else
      error('invalid preconditioner!')
   end
%
   % zero initial guess
   x0=zeros(size(fgal));
   tic %%start timing
   if itmeth==1, %PCG
      fprintf('\nPCG iteration ...\n');
      if precon<=2, 
         [x_it,flag,relres,iter,resvec] = pcg(Agal,fgal,tol,maxit,M1,M2,x0);
      elseif  precon==3
          [x_amg,flag,relres,iter,resvec] = ...
             pcg(Agal,fgal,tol,maxit, @amg_v_cycle, [], x0, amg_grid, amg_smoother); 
      end
   elseif itmeth==2, %MINRES
      fprintf('\nMINRES iteration ...\n');
      if precon<=2, 
         [x_it,flag,relres,iter,resvec] = minres(Agal,fgal,tol,maxit,M1,M2,x0);
      elseif  precon==3
           [x_it,flag,relres,iter,resvec] = ...
              minres(Agal,fgal,tol,maxit,@amg_v_cycle, [], x0, amg_grid, amg_smoother);
      end
   else
      error('invalid iterative method!')
   end
   etoc = toc;
%
else, error('Oops. solvers for this PDE are not available yet')
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print and plot results
if flag ==0,
   % successful convergence
   fprintf('convergence in %3i iterations\n',iter)
   nr0=resvec(1);
   fprintf('\n    k  log10(||r_k||/||r_0||)   \n')
   for its=1:iter+1,
      fprintf('%5i %16.4f \n', its-1, log10(resvec(its)/nr0));
   end
   fprintf('Bingo!\n')
   fprintf('\n  %9.4e seconds\n\n\n',etoc)  
   %%% plot residuals
   resplot(resvec)
else
   nr0=resvec(1);
   fprintf('\n    k  log10(||r_k||/||r_0||)   \n')
   for its=1:iter+1,
      fprintf('%5i %16.4f \n', its-1, log10(resvec(its)/nr0));
   end
   fprintf('iteration aborted! Iteration returned with flag equal to  %2i \n',flag)
   %%% plot residuals
   resplot(resvec)
end