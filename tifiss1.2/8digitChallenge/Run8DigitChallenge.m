%RUN8DIGITCHALLENGE  posed by L. N. Trefethen on NA-digest (Windows version)
%   TIFISS scriptfile: AB; 08 January 2019
% Copyright (c) 2018 D.J. Silvester, A. Bespalov

  close all;
  gohome;
  
  fprintf('\nRun 8 digit challenge ...');
  fprintf('\nL-shaped domain, harmonic solution, Dirichlet b.c.');
  fprintf('\n');
   system('copy .\diffusion\test_problems\unit_adiff.m .\diffusion\specific_adiff.m');
   system('copy .\diffusion\test_problems\zeros_gradcoeff.m .\diffusion\specific_gradcoeff.m');
   system('copy .\diffusion\test_problems\zero_rhs.m .\diffusion\specific_rhs.m');
   system('copy .\8digitChallenge\bc_ex8.m .\diffusion\specific_bc.m');

   dom_type = 2; sn=8;

% Call the main driver
  compute8digits;
nit=length(P2valpoint);
pointestimate=[[1:nit];P2valpoint];
fprintf('\nBingo!')
fprintf('\n  itn   point value\n')
fprintf('   %2i %15.11f \n',pointestimate)
fprintf('\n')

% end scriptfile  
