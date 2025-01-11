%HELPME_DIFF   anisotropic diffusion problem interactive help
%   TIFISS scriptfile: AB; 11 February 2019.
% Copyright (c) 2016 D.J. Silvester, Alex Bespalov, Qifeng Liao

fprintf(' \n');
fprintf(' To generate and solve a diffusion equation with variable coefficient\n');
fprintf(' over a square or an L-shaped domain, run the driver: <strong>diff_testproblem</strong>\n');
fprintf(' \n');
fprintf(' Typing <CR> when prompted for input automatically gives the default choice.\n');
fprintf(' \n');
fprintf(' Dirichlet boundary conditions and forcing terms are set in the user-defined functions:\n');
fprintf('     <strong>/diffusion/specific_bc.m</strong>  and  <strong>/diffusion/specific_rhs.m</strong>\n');
fprintf(' Diffusion tensor (diagonal) coefficients are set in the user-defined functions:\n');
fprintf('     <strong>/diffusion/specific_adiff.m</strong>  and  <strong>/diffusion/specific_gradcoeff.m</strong>');
fprintf(' \n\n');
fprintf(' One can also generate a self-adaptive solution by running:\n'); 
fprintf(' <strong>adapt_diff_testproblem</strong>  or  <strong>goafem_testproblem</strong> (goal-oriented)');
fprintf(' \n\n');
fprintf(' To generate and solve a collection of reference diffusion problems over\n');
fprintf(' the square, L-shaped, circle-shaped, and punched ticket domains,\n');
fprintf(' run the driver: <strong>diff_refproblem</strong>');
fprintf(' \n');


