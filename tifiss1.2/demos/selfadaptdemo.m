%SELFADAPTDEMO adaptive refinement
%   TIFISS scriptfile: DJS; 5 March 2017.
%Copyright (c) 2017 by D.J. Silvester, Alex Bespalov
tifiss, pause(1)
close all
fprintf('1. square domain | isotropic diffusion ... \n'), pause(3)
batchmode('adaptP1'),
fprintf('\nCHECK OUT the computed solution \n'), pause(5)
fprintf('CHECK OUT the final grid \n'), figure(2), pause(5)
fprintf('CHECK OUT error reduction path \n'), figure(3)
fprintf(' (Type any character to continue.)\n'), pause
fprintf('2. square domain | strong anisotropy ... \n'), pause(3)
batchmode('adaptP2'),
fprintf('\nCHECK OUT the computed solution \n'), pause(5)
fprintf('CHECK OUT the final grid \n'), figure(2)
fprintf(' (Type any character to continue.)\n'), pause
fprintf('5. L-shaped domain | refinement of singularity ...\n'), pause(3)
batchmode('adaptP5'),
fprintf('\nCHECK OUT the computed solution \n'), pause(5)
fprintf('CHECK OUT the final grid \n'), figure(2), pause(5)
fprintf('CHECK OUT error reduction path \n'), figure(3)
fprintf('Voila! end of diffusion demo. \n');
