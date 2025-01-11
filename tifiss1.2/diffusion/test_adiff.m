%TEST_ADIFF   solve (an)isotropic problem in half L-shaped domain
%   TIFISS scriptfile: DJS; 5 March 2017.
% Copyright (c) 2017 A. Bespalov, D.J. Silvester

%% define geometry
pde=1; domain=12;
test_domain
load test_grid
figure(10)
imeshplot(xy,evt,bound,eboundt);
% validate mesh
[eex,tve,els] = tedgegen(xy,evt);
%
%-------------------P1 approximation only
%% set up matrices
qmethod=1;
[A,M,f] = femp1_adiff(xy,evt);

%% boundary conditions
[Agal,fgal] = nonzerobc(A,f,xy,bound);

%% compute solution
tic
fprintf('solving linear system ...  ')
x_gal = Agal\fgal;
fprintf('done\n')
etoc=toc; fprintf('Galerkin system solved in  %8.3e seconds\n',etoc) 

%% Test alternative a posteriori error estimate(s)
pestim = default('Error estimation: linear/quadratic bubble functions 1/2? (default 1)',1);
      
if pestim == 1
       
% LINEAR MIDPOINT HAT FUNCTIONS
      fprintf('\n<strong>diffpost using 3 edge midpoint linear hat functions...</strong>\n'); 
      tic; [elerr_p,fe,ae] = diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,x_gal); toc;
           %fprintf('sanity check %10.4e\n',norm(elerr_p-elerr_x,inf))
      [err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
      
elseif pestim == 2
            
% QUADRATIC MIDPOINT BUBBLE FUNCTIONS
      fprintf('\n<strong>diffpost using 4 quadratic bubble functions...</strong>\n');
      tic; [elerr_p,fe,ae] = diffpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,x_gal); toc;    
          %fprintf('sanity check %10.4e\n',norm(elerr_p-elerr_x,inf))
      [err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
      
else
       error('Value not valid! Choose between either 1 or 2!');
end
     
%  terrplot(qmethod,x_gal,elerr_p,evt,xy,x,y,101,[xmin,xmax,ymin,ymax]);
fprintf('<strong>Estimated energy error: %10.4e</strong>\n',norm(elerr_p,2));
