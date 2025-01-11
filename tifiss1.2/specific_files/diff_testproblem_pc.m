%diff_testproblem_pc   sets up test examples  (Windows version)
%   TIFISS scriptfile: DJS; 21 December 2017.
%Copyright (c) 2016 by D.J. Silvester, Alex Bespalov, Qifeng Liao

gohome
clear variables
close all
fprintf('\nspecification of reference diffusion problem.')
fprintf('\nchoose specific example');
fprintf('\n     1  Square domain, unit coefficients')
fprintf('\n     2  Square domain, anisotropic coefficients')
fprintf('\n     3  Square domain, variable coefficients')
fprintf('\n     4  L-shaped domain, constant load function')
fprintf('\n     5  L-shaped domain, variable coefficients, analytic solution')
fprintf('\n')

sn = default('',1);

if sn==1
   system('copy .\diffusion\test_problems\unit_rhs.m .\diffusion\specific_rhs.m');
   system('copy .\diffusion\test_problems\zero_bc.m .\diffusion\specific_bc.m');
   system('copy .\diffusion\test_problems\unit_adiff.m .\diffusion\specific_adiff.m');
   system('copy .\diffusion\test_problems\zeros_gradcoeff.m .\diffusion\specific_gradcoeff.m');
   square_adiff
elseif sn==2
   system('copy .\diffusion\test_problems\unit_rhs.m .\diffusion\specific_rhs.m');
   system('copy .\diffusion\test_problems\zero_bc.m .\diffusion\specific_bc.m');
   system('copy .\diffusion\test_problems\anisotropic_adiff.m .\diffusion\specific_adiff.m');
   system('copy .\diffusion\test_problems\zeros_gradcoeff.m .\diffusion\specific_gradcoeff.m');
   square_adiff
elseif sn==3
    system('copy .\diffusion\test_problems\rhs_ex6.m .\diffusion\specific_rhs.m');
    system('copy .\diffusion\test_problems\zero_bc.m .\diffusion\specific_bc.m');
    system('copy .\diffusion\test_problems\adiff_ex6.m .\diffusion\specific_adiff.m');
    system('copy .\diffusion\test_problems\gradcoeff_ex6.m .\diffusion\specific_gradcoeff.m');
    square_adiff
elseif sn==4
   system('copy .\diffusion\test_problems\unit_rhs.m .\diffusion\specific_rhs.m');
   system('copy .\diffusion\test_problems\zero_bc.m .\diffusion\specific_bc.m');
   system('copy .\diffusion\test_problems\unit_adiff.m .\diffusion\specific_adiff.m');
   system('copy .\diffusion\test_problems\zeros_gradcoeff.m .\diffusion\specific_gradcoeff.m');
   ell_adiff
elseif sn==5
   system('copy .\diffusion\test_problems\rhs_ex7.m .\diffusion\specific_rhs.m');
   system('copy .\diffusion\test_problems\bc_ex7.m .\diffusion\specific_bc.m');
   system('copy .\diffusion\test_problems\adiff_ex7.m .\diffusion\specific_adiff.m');
   system('copy .\diffusion\test_problems\gradcoeff_ex7.m .\diffusion\specific_gradcoeff.m');
   ell_adiff
else
   error('Oops.. reference problem datafile not found!')
end  % end if
