function sm_swi_jobs(reg, res)

% sm_swi_jobs.m
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if strcmp('sir',res)
    [lmask, ~] = loadsir(['landmasks/' reg '.sir.lmask']);
    [num_rows, ~] = size(lmask);    
elseif strcmp('grd',res)
    [lmask, ~] = loadsir(['landmasks/' reg '.grd.lmask']);
    [num_rows, ~] = size(lmask);        
else
    error('Resolution option not supported!');
end
for ii = 1:num_rows %[77 78 79 80 81 82 83 84 363 537 538 539 540 541 542 543 814 ]%

    PATH = getenv('PATH');
    LD_LIBRARY_PATH = getenv('LD_LIBRARY_PATH');
    cd('/home/lindell/research/soil_moisture');

    filename = ['slurm_jobs/' 'swi_' reg '_' sprintf('%04d',ii) '.job'];

    filetext = ...
    [sprintf('%s\n','#!/usr/bin/env bash') ...
    sprintf('%s\n',['#SBATCH -J' ' swi_' reg '_' sprintf('%04d',ii)  '.job' ' --time=1:00:00 --mem-per-cpu=2048 --ntasks=1 --output=/home/lindell/research/soil_moisture/slurm_out/slurm-%%j.out --mail-user=lindell@mers.byu.edu --mail-type=FAIL']) ...
    sprintf('%s\n','hostname') ...
    sprintf('%s\n',['export PATH=', PATH]) ...
    sprintf('%s\n',['export LD_LIBRARY_PATH=', LD_LIBRARY_PATH]) ...
    sprintf('%s\n','cd /home/lindell/research/soil_moisture') ...
    sprintf('%s\n',['matlab -nojvm -nodisplay -r "cd(''/home/lindell/research/soil_moisture''),addpath(''/auto/users/long/matcode:/home/lindell/Documents/MATLAB/Apps:/home/lindell/Documents/MATLAB''),sm_gen_swi(' '''' reg '''' ',' num2str(ii) ',' '''' res '''' '),quit()"'])];


    fid = fopen(filename,'w');
    fprintf(fid,filetext);
    fclose(fid);

    command = ['sbatch ' filename];
    system(command);
    
end