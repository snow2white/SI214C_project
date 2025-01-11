%DIFFPROBLEMDEMO diffusion test problems
%   TIFISS scriptfile: DJS; 12 February 2016.
%Copyright (c) 2016 by D.J. Silvester, Qifeng Liao
tifiss, pause(1)
close all
fprintf('\n1. square domain | unit coefficients ... \n'), pause(3)
batchmode('P1'),
fprintf('\nCHECK OUT the error in the computed solution\n'),
fprintf(' (Type any character to continue.)\n'), pause
fprintf('2. square domain | strong anisotropy ... \n'), pause(3)
batchmode('P2'),
fprintf('\nCHECK OUT the computed solution \n')
fprintf(' (Type any character to continue.)\n'), pause
fprintf('3.circular domain | exact quadratic solution ... \n'), pause(3)
batchmode('P3'),
fprintf('\nCHECK OUT the perfectly computed solution\n'),
fprintf(' (Type any character to continue.)\n'), pause
fprintf('4. obstacle domain | smooth solution\n'), pause(3)
fprintf('be patient .. the grid generation takes a while ... \n')
batchmode('P4'),
fprintf('\nCHECK OUT the computed solution\n'),
fprintf(' (Type any character to continue.)\n'), pause
fprintf('5.L-shaped domain | refinement of singularity ...\n'), pause(3)
batchmode('P5'),
fprintf('\nCHECK OUT the computed solution \n'),
fprintf(' (Type any character to continue.)\n'), pause
% iterative solver
fprintf('finally CHECK the AMG solver convergence ...\n'), pause(3)
batchmode('itsolve_amg')
fprintf('Voila! end of diffusion demo. \n'), figure(19)
