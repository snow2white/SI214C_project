%INSTALL_PC_X sets up TIFISS on Windows computer
%    TIFISS scriptfile: LR; 24 January 2019, DJS; 24 April 2020
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

if ispc %strncmp(computer,'PCWIN',3)
   gohome;
   fprintf('\nInstalling PC-OS specific files.\n');
%
   delete .\diffusion\test_problems\diff_testproblem.m;
   system('copy .\specific_files\diff_testproblem_pc.m          .\diffusion\test_problems\diff_testproblem.m');
   delete .\diffusion\test_problems\diff_refproblem.m;
   system('copy .\specific_files\diff_refproblem_pc.m           .\diffusion\test_problems\diff_refproblem.m');
   delete .\diffusion\test_problems\adapt_diff_testproblem.m;
   system('copy .\specific_files\adapt_diff_testproblem_pc.m    .\diffusion\test_problems\adapt_diff_testproblem.m');
%
% GOAFEM
   delete .\goafem\test_problems\goafem_testproblem.m;
   system('copy .\specific_files\goafem_testproblem_pc.m  .\goafem\test_problems\goafem_testproblem.m');
% 8digits
   delete .\8digitChallenge\Run8DigitChallenge.m'
   system('copy .\specific_files\Run8DigitChallenge_pc.m      .\8digitChallenge\Run8DigitChallenge.m');
%
   %
   %
   if exist('stoch_diffusion','dir')
       gohome; cd stoch_diffusion;
       delete .\test_problems\stoch_diff_testproblem.m;
       system('copy .\test_problems\stoch_diff_testproblem_pc.m   .\test_problems\stoch_diff_testproblem.m');
       delete .\test_problems\stoch_adapt_testproblem.m;
       system('copy .\test_problems\stoch_adapt_testproblem_pc.m  .\test_problems\stoch_adapt_testproblem.m');
       % STOCH_GOAFEM
       delete .\test_problems\stoch_goafem_testproblem.m;
       system('copy .\test_problems\stoch_goafem_testproblem_pc.m  .\test_problems\stoch_goafem_testproblem.m');
   end
   %
   fprintf('Done.\n');
else
   fprintf('\nOops.. run install_pc to get PC_OS files on a Unix machine\n');
end
% end scriptfile
