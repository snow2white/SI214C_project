%GOAFEM_TEST_PROBLEM sets up goal-oriented adaptive FEM examples (Windows version)
%   TIFISS scriptfile: MR; 04 July 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  gohome;

  fprintf('\nSpecification of goal-oriented example for adaptive FEM')
  fprintf('\nChoose specific example:')
  fprintf('\n   1.  Square domain (0,1)^2, symmetric Mommer-Stevenson example;');
  fprintf('\n   2.  L-shaped domain, symmetric Mommer-Stevenson example;');
  fprintf('\n   3.  L-shaped domain, nonsymmetric Mommer-Stevenson example;');
  fprintf('\n   4.  Crack domain, pointwise estimation;\n');
    
% Running by default the Mommer-Stevenson-type example
  sn = default('',1);
  
% NOTE:
% sn: example problem
%
% dom_type: domain type
% - 1 unit square domain (0,1)^2
% - 2 L-shaped domain    (-1,1)^2 \ (-1,0]^2
% - 3 large crack domain (-1,1)^2 \ (-1,0)x{0} (crack on the left)

% Set up data for specific problems
  if sn==1
      % Mommer-Stevenson-type example on square domain
      system('copy .\goafem\test_problems\goafem_unit_coeff.m       .\goafem\goafem_specific_coeff.m');
      system('copy .\goafem\test_problems\goafem_zero_gradcoeff.m   .\goafem\goafem_specific_gradcoeff.m');
      system('copy .\goafem\test_problems\goafem_zero_bc.m          .\goafem\goafem_specific_bc.m');
      system('copy .\goafem\test_problems\goafem_zero_L2rhs.m       .\goafem\goafem_specific_L2rhs.m');
      system('copy .\goafem\test_problems\goafem_zero_L2goal.m      .\goafem\goafem_specific_L2goal.m');
      system('copy .\goafem\test_problems\goafem_ms09_H1rhs.m       .\goafem\goafem_specific_H1rhs.m');
      system('copy .\goafem\test_problems\goafem_ms09_H1goal.m      .\goafem\goafem_specific_H1goal.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1rhs.m    .\goafem\goafem_specific_divH1rhs.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1goal.m   .\goafem\goafem_specific_divH1goal.m');
      dom_type = 1;
  
  elseif sn==2     
      % Symmetric Mommer-Stevenson-type example on L-shaped domain
      system('copy .\goafem\test_problems\goafem_unit_coeff.m       .\goafem\goafem_specific_coeff.m');
      system('copy .\goafem\test_problems\goafem_zero_gradcoeff.m   .\goafem\goafem_specific_gradcoeff.m');
      system('copy .\goafem\test_problems\goafem_zero_bc.m          .\goafem\goafem_specific_bc.m');
      system('copy .\goafem\test_problems\goafem_unit_L2rhs.m       .\goafem\goafem_specific_L2rhs.m');
      system('copy .\goafem\test_problems\goafem_unit_L2goal.m      .\goafem\goafem_specific_L2goal.m');
      system('copy .\goafem\test_problems\goafem_mod_ms09_H1rhs.m   .\goafem\goafem_specific_H1rhs.m');
      system('copy .\goafem\test_problems\goafem_mod_ms09_H1goal.m  .\goafem\goafem_specific_H1goal.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1rhs.m    .\goafem\goafem_specific_divH1rhs.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1goal.m   .\goafem\goafem_specific_divH1goal.m');
      dom_type = 2;
 
  elseif sn==3
      % Nonsymmetric Mommer-Stevenson-type example on L-shaped domain
      system('copy .\goafem\test_problems\goafem_unit_coeff.m       .\goafem\goafem_specific_coeff.m');
      system('copy .\goafem\test_problems\goafem_zero_gradcoeff.m   .\goafem\goafem_specific_gradcoeff.m');
      system('copy .\goafem\test_problems\goafem_zero_bc.m          .\goafem\goafem_specific_bc.m');
      system('copy .\goafem\test_problems\goafem_unit_L2rhs.m       .\goafem\goafem_specific_L2rhs.m');
      system('copy .\goafem\test_problems\goafem_zero_L2goal.m      .\goafem\goafem_specific_L2goal.m');
      system('copy .\goafem\test_problems\goafem_zero_H1rhs.m       .\goafem\goafem_specific_H1rhs.m');
      system('copy .\goafem\test_problems\goafem_mod_ms09_H1goal.m  .\goafem\goafem_specific_H1goal.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1rhs.m    .\goafem\goafem_specific_divH1rhs.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1goal.m   .\goafem\goafem_specific_divH1goal.m');      
      dom_type = 2;

  elseif sn==4 
      % Pointwise estimation on large crack-domain (-1,1)^2 \ (-1,0)x{0}
      system('copy .\goafem\test_problems\goafem_unit_coeff.m       .\goafem\goafem_specific_coeff.m');
      system('copy .\goafem\test_problems\goafem_zero_gradcoeff.m   .\goafem\goafem_specific_gradcoeff.m');
      system('copy .\goafem\test_problems\goafem_zero_bc.m          .\goafem\goafem_specific_bc.m');
      system('copy .\goafem\test_problems\goafem_unit_L2rhs.m       .\goafem\goafem_specific_L2rhs.m');
      system('copy .\goafem\test_problems\goafem_po09_L2goal.m      .\goafem\goafem_specific_L2goal.m');
      system('copy .\goafem\test_problems\goafem_zero_H1rhs.m       .\goafem\goafem_specific_H1rhs.m');
      system('copy .\goafem\test_problems\goafem_zero_H1goal.m      .\goafem\goafem_specific_H1goal.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1rhs.m    .\goafem\goafem_specific_divH1rhs.m');
      system('copy .\goafem\test_problems\goafem_zero_divH1goal.m   .\goafem\goafem_specific_divH1goal.m');
      dom_type = 3;
  else
      error('Oops... reference problem datafile not found!');
  end % end if

% Calling the main driver
  goafem_diffusion_main;
  
% end scriptfile