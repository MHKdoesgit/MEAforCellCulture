

function InesRunKilosort()

ksfunctionpath = 'C:\MEA\Matlab_functions\KiloSortMEA';

disp('!!!!!MAY THE FORCE BE WITH YOU!!!!!');

addpath(genpath (ksfunctionpath));

batch_run_ks;

rmpath(ksfunctionpath);

disp('SORTING IS FINITO');


end

