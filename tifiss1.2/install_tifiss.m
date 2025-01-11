%INSTALL_TIFISS installation instructions
%   TIFISS scriptfile: DS; 20 April 2020
%Copyright (c) 2016 by D.J. Silvester, Qifeng Liao
%
fprintf(' To install the toolbox, the script-file\n');
fprintf(' gohome.m must be edited to reflect the correct path\n');
fprintf(' to the TIFISS home directory on the installed computer: \n')
fprintf(' \n %s\n\n',pwd);
%currentdir=pwd; fprintf(['\n ',currentdir,' \n\n']);
fprintf(' This step only needs to be done once.\n');
fprintf(' (Type any character to continue.)')
pause;


fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
fprintf(' \n');
fprintf(' Once T-IFISS is installed, for all subsequent uses the MATLAB search path\n');
fprintf(' must include the TIFISS subdirectories. This can be done by running \n');
fprintf(' the script-file setpath.m from the TIFISS home directory.\n');
fprintf(' (Type any character to continue.)')
pause;

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
fprintf(' \n');
fprintf(' To test the installation, you can regenerate the saved sessions \n')
fprintf(' in the included diary files; for example, type \n')
fprintf(' "edit testinstall.txt", "edit testadapt_ex5.txt" or "edit goafem_testadapt_ex1.txt"\n')

