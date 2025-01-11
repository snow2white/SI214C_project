%SETPATH sets TIFISS search path
%   TIFISS scriptfile: DJS; 12 January 2016.
%  Copyright (c) 2016 D.J. Silvester, Qifeng Liao


% initilise matlab path for ./tifiss/ directories
gohome, 
warning('off', 'all')
addpath(genpath(pwd),'-end')
cd matlab704, rmpath(pwd), gohome
cd matlabpre74, rmpath(pwd), gohome
%fprintf('Resolving these path conflicts ... done')
warning('on','all')
%

% fix iterapp bug (affects pcg & minres)
if strncmp(version,'7.0.4',5), 
cd matlab704, addpath(pwd), gohome
end

% include old version of ilu 
if strncmp(version,'6.5',3) | strncmp(version,'7.0',3) | strncmp(version,'7.1.',4) | ...
    strncmp(version,'7.2',3) | strncmp(version,'7.3',3), 
cd matlabpre74, addpath(pwd), gohome
end

% fix MATLAB UMFPACK bug
if strncmp(version,'6.5',3) | strncmp(version,'7.0',3), 
spparms('default'),spparms('piv_tol',1), end

