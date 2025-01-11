%HELPME   TIFISS interactive help facility
%   TIFISS scriptfile: DJS; 5 March 2017.
%Copyright (c) 2017 by D.J. Silvester, Alex Bespalov, Qifeng Liao, Leonardo Rocchi

fprintf(' \n');
fprintf(' T-IFISS\n')
fprintf('Copyright (c) 2017 by D.J. Silvester and Alex Bespalov and Qifeng Liao and Leonardo Rocchi\n')
%fprintf(' \n');
fprintf('To install the toolbox, run the script-file install_tifiss.m\n');
%fprintf(' \n');

fprintf(' (Type any character to continue.)')
pause;

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
fprintf('\nTry running the software demo(s) to get started\n'); help diffproblemdemo, help selfadaptdemo

fprintf('\n For further information on\n');
fprintf('    new features                                      enter 0\n');
fprintf('    solving a deterministic diffusion problem               1\n');
fprintf('    solving a stochastic diffusion problem                  2\n');
fprintf('    exploring preconditioned Krylov subspace solvers        99\n\n');
hlp=default('    Help topic',-1);
if hlp==1, gohome; cd diffusion; helpme_diff;
elseif hlp==2,
if exist('stoch_diffusion')
 gohome; cd stoch_diffusion; helpme_stoch_diffusion;
else, fprintf('\n Oops .. stochastic directory is not included in vanilla T-IFISS\n'); end
elseif hlp==99, gohome; cd solvers; helpme_it;
elseif hlp==0,
gohome; type release.txt
end
gohome;
fprintf(' \n');

