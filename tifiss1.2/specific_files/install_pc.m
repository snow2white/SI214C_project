%INSTALL_PC sets up Windows TIFISS on UNIX computer
%    TIFISS scriptfile: DJS 18 December 2018; 18 August 2019
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao

if isunix %strncmp(computer,'MACI',4)
    gohome;
    fprintf('\nInstalling PC-OS specific files.\n')
%
system('/bin/rm ./diffusion/test_problems/diff_testproblem.m');
system('/bin/cp ./specific_files/diff_testproblem_pc.m ./diffusion/test_problems/diff_testproblem.m');
system('/bin/rm ./diffusion/test_problems/diff_refproblem.m');
system('/bin/cp ./specific_files/diff_refproblem_pc.m ./diffusion/test_problems/diff_refproblem.m');
system('/bin/rm ./diffusion/test_problems/adapt_diff_testproblem.m');
system('/bin/cp ./specific_files/adapt_diff_testproblem_pc.m ./diffusion/test_problems/adapt_diff_testproblem.m');
%
% GOAFEM
system('/bin/rm ./goafem/test_problems/goafem_testproblem.m');
system('/bin/cp ./specific_files/goafem_testproblem_pc.m  ./goafem/test_problems/goafem_testproblem.m');
%
% 8digits
system('/bin/rm ./8digitChallenge/Run8DigitChallenge.m');
system('/bin/cp ./specific_files/Run8DigitChallenge_pc.m ./8digitChallenge/Run8DigitChallenge.m');
%
if exist('stoch_diffusion')
gohome; cd stoch_diffusion;
system('/bin/rm ./test_problems/stoch_diff_testproblem.m');
system('/bin/cp ./test_problems/stoch_diff_testproblem_pc.m  ./test_problems/stoch_diff_testproblem.m');
system('/bin/rm ./test_problems/stoch_adapt_testproblem.m');
system('/bin/cp ./test_problems/stoch_adapt_testproblem_pc.m  ./test_problems/stoch_adapt_testproblem.m');
system('/bin/rm  ./test_problems/stoch_goafem_testproblem.m');
system('/bin/cp ./test_problems/stoch_goafem_testproblem_pc.m  ./test_problems/stoch_goafem_testproblem.m');
end
%
   fprintf('Done.\n')
end

% end scriptfile
