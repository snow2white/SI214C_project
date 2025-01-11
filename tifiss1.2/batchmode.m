function batchmode(testproblem)
%BATCHMODE enables batch processing for TIFISS testproblem
%   batchmode(testproblem);
%   input
%          testproblem  character string naming the testproblem
%                       must have the form "*_batch".m where "*" begins with
%                       "P"       for diffusion problems
%                      "adapt"    for adptive refinement run
%   side effect
%          If batchmode terminates prematurely because of an error or
%          execution of cntl-C, interactive input with IFISS may not
%          work correctly.  This is fixed by typing "activemode".
%
%
%   IFISS function: HCE, DJS; 6 March 2017.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 

global BATCH FID

% file containing input data 
batchfile=[testproblem,'_batch.m'];
[FID,message]=fopen(batchfile,'r');
% Error return on nonexistent or misnamed input file
if strcmp(message,'')~=1
   error(['INPUT FILE ERROR: ' message])
else
   disp(['Working in batch mode from data file ' batchfile])
end
if ~(strncmp(testproblem,'P',1) |  strncmp(testproblem,'itsolve',7)| strncmp(testproblem,'adapt',5)),
    errmsg = 'INPUT FILE ERROR:\n';
    errmsg = [errmsg,'   Batch input filenames must have the form "*_batch.m"'];
    errmsg = [errmsg,' where "*" begins with\n'];
    errmsg = [errmsg,'   "P" for generation of diffusion problems\n'];
    errmsg = [errmsg,'   "adapt" for adaptive refinement run.'];
    errmsg = [errmsg,'   "itsolve" for iterative solution of linear systems.']
    error('BATCH:err',errmsg);    
end  

% batch run
% switch to activate batch mode (off/on 0/1) (see "default.m")
BATCH=1;

% run appropriate driver
   if strncmp(testproblem,'P',1)  
      diff_refproblem
   gohome, cd datafiles, save batchrun.mat
   elseif strncmp(testproblem,'adapt',5)
      adapt_diff_testproblem
   gohome, cd datafiles, save adaptrun.mat
   elseif strncmp(testproblem,'itsolve',7)
   load batchrun.mat
   it_solve
   gohome, cd datafiles
   save batchrun_itsolve.mat
   end


% switch back to interactive mode
activemode
return
