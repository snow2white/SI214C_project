%ADAPT_DIFF_TESTPROBLEM sets up adaptive mesh test examples (Unix version)
%   TIFISS scriptfile: LR; 22 June 2018, DJS; 16 August 2019
%Copyright (c) 2016 by D.J. Silvester, A. Bespalov

  clear variables, close all;
  gohome;
  
  fprintf('\nSpecification of reference diffusion problem.');
  fprintf('\nChoose specific example');
  fprintf('\n   1.  Square domain,   unit coefficients');
  fprintf('\n   2.  Square domain,   anisotropic coefficients');
  fprintf('\n   3.  Square domain,   variable coefficients');
  fprintf('\n   4.  L-shaped domain, constant load function');
  fprintf('\n   5.  L-shaped domain, variable coefficients, analytic solution');
  fprintf('\n   6.  Crack domain,    variable coefficients');
  fprintf('\n');

% Running example 1 by default
  sn = default('',1);

% NOTE:
% sn: example problem
%
% dom_type: domain type
% - 1 unit square domain (0,1)^2
% - 2 L-shaped domain    (-1,1)^2 \ (-1,0]^2
% - 3 large crack domain (-1,1)^2 \ (-1,0)x{0} (crack on the left)

% Example problem
  if sn==1
      system('/bin/cp ./diffusion/test_problems/unit_adiff.m         ./diffusion/specific_adiff.m');
      system('/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m    ./diffusion/specific_gradcoeff.m');
      system('/bin/cp ./diffusion/test_problems/unit_rhs.m           ./diffusion/specific_rhs.m');
      system('/bin/cp ./diffusion/test_problems/zero_bc.m            ./diffusion/specific_bc.m');
      dom_type = 1;
  elseif sn==2
      system('/bin/cp ./diffusion/test_problems/anisotropic_adiff.m  ./diffusion/specific_adiff.m');
      system('/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m    ./diffusion/specific_gradcoeff.m');
      system('/bin/cp ./diffusion/test_problems/unit_rhs.m           ./diffusion/specific_rhs.m');
      system('/bin/cp ./diffusion/test_problems/zero_bc.m            ./diffusion/specific_bc.m');
      dom_type = 1;
  elseif sn==3
      system('/bin/cp ./diffusion/test_problems/adiff_ex6.m          ./diffusion/specific_adiff.m');
      system('/bin/cp ./diffusion/test_problems/gradcoeff_ex6.m      ./diffusion/specific_gradcoeff.m');
      system('/bin/cp ./diffusion/test_problems/rhs_ex6.m            ./diffusion/specific_rhs.m');
      system('/bin/cp ./diffusion/test_problems/zero_bc.m            ./diffusion/specific_bc.m');
      dom_type = 1;
  elseif sn==4
      system('/bin/cp ./diffusion/test_problems/unit_adiff.m         ./diffusion/specific_adiff.m');
      system('/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m    ./diffusion/specific_gradcoeff.m');
      system('/bin/cp ./diffusion/test_problems/unit_rhs.m           ./diffusion/specific_rhs.m');
      system('/bin/cp ./diffusion/test_problems/zero_bc.m            ./diffusion/specific_bc.m');
      dom_type = 2;
  elseif sn==5    
      % Exact singular solution on L-shape with non-homogeneous bcs:
      % Solution: u(r,theta) = r^(2/3) * sin( (2*theta + pi) / 3 )
      % This solution solves the following two equivalent problems:
      % 1. div( (1 - r^2/4) * \grad(u(r,theta)) ) = u(r,theta) / 3
      % 2. Lap( u(r,theta) ) = 0
      % Choose the problem (1/2)?
      problem = 1;
      if problem == 1
      system('/bin/cp ./diffusion/test_problems/adiff_ex7.m          ./diffusion/specific_adiff.m');
      system('/bin/cp ./diffusion/test_problems/gradcoeff_ex7.m      ./diffusion/specific_gradcoeff.m');
      system('/bin/cp ./diffusion/test_problems/rhs_ex7.m            ./diffusion/specific_rhs.m');
      system('/bin/cp ./diffusion/test_problems/bc_ex7.m             ./diffusion/specific_bc.m');
      else%if problem == 2
      system('/bin/cp ./diffusion/test_problems/unit_adiff.m         ./diffusion/specific_adiff.m');
      system('/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m    ./diffusion/specific_gradcoeff.m');
      system('/bin/cp ./diffusion/test_problems/zero_rhs.m           ./diffusion/specific_rhs.m');
      system('/bin/cp ./diffusion/test_problems/bc_ex7.m             ./diffusion/specific_bc.m');
      end
      dom_type = 2;
  elseif sn==6
      system('/bin/cp ./diffusion/test_problems/unit_adiff.m         ./diffusion/specific_adiff.m');
      system('/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m    ./diffusion/specific_gradcoeff.m');
      system('/bin/cp ./diffusion/test_problems/unit_rhs.m           ./diffusion/specific_rhs.m');
      system('/bin/cp ./diffusion/test_problems/zero_bc.m            ./diffusion/specific_bc.m');
      dom_type = 3;
  else
      error('Oops.. reference problem datafile not found!');
  end  % end if

% Call the main driver
  adiff_adaptive_main;

% end scriptfile
