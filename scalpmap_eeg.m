% This file is for eeg-data classification
% All MATLAB functions are available in BBCI toolbox
% Some minor code modifications might be applied
% We do not guarantee all of functions works properly in your platform
% For more tutorials, visit BBCI toolbox (https://github.com/bbci/bbci_public)

% specify your eeg data directory (EegMyDataDir) and temporary directory (TemDir)
startup_bbci_toolbox('DataDir',EegMyDataDir,'TmpDir','/tmp/');
BTB.History = 0; % to aviod error for merging cnt
fs = 200 % downsampling rate: 20 Hz

%% initial parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subdir_list = {'subject1'}
basename_list = {'imagery1','arithmetic1','imagery1','arithmetic2','imagery1','arithmetic3'};
stimDef.eeg= {16, 32;'condition1','condition2'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify discriminative frequency band. Use alpha band for example.
band_csp.imag = [8 14]; % alpha-band for motor imagery

% motor imagery
loadDir = fullfile(EegMyDataDir,subdir_list{1});
cd(loadDir);
load cnt; load mrk; load mnt;
cd(WorkingDir);

cnt_temp = cnt; mrk_temp = mrk; % backup
clear cnt mrk;

[cnt.imag, mrk.imag] = proc_appendCnt({cnt_temp{1}, cnt_temp{3}, cnt_temp{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
[b,a]= butter(3, band_csp.imag/fs*2);
cnt.imag = proc_filt(cnt.imag, b, a);

% Specify discriminative frequency band. Use alpha band for example.
band_csp.ment = [8 14]; % alpha-band for motor mentery

% motor mentery
loadDir = fullfile(EegMyDataDir,subdir_list{1});
cd(loadDir);
load cnt; load mrk; load mnt;
cd(WorkingDir);

cnt_temp = cnt; mrk_temp = mrk; % backup
clear cnt mrk;

[cnt.ment, mrk.ment] = proc_appendCnt({cnt_temp{2}, cnt_temp{4}, cnt_temp{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}}); % merged mental arithmetic cnts
[b,a]= butter(3, band_csp.ment/fs*2);
cnt.ment = proc_filt(cnt.ment, b, a);

ival_epo = [-10 25]*1000;
ival_scalps = [-8 -5; -3 0; 2 5; 7 10; 12 15; 17 20]*1000;
ival_base = [-3 0]*1000;

% segmentation
% motor imagery
epo.imag = proc_segmentation(cnt.imag, mrk.imag, ival_epo);
epo.className = {'LMI','RMI'};
epo_car.imag = proc_commonAverageReference(epo.imag);
epo_car.imag = proc_envelope(erd_car.imag, 'MovAvgMsec', 500);
epo_car_r.imag = proc_rSquareSigned(erd_car.imag);

fig_set(1)
plot_scalpEvolution(erd_car_r, mnt_grand, ival_scalps, defopt_scalp_r);

% mental arithmetic
epo.ment = proc_segmentation(cnt.ment, mrk.ment, ival_epo);
epo.className = {'LMI','RMI'};
epo_car.ment = proc_commonAverageReference(epo.ment);
epo_car.ment = proc_envelope(erd_car.ment, 'MovAvgMsec', 500);
epo_car_r.ment = proc_rSquareSigned(erd_car.ment);

fig_set(2)
plot_scalpEvolution(erd_car_r, mnt_grand, ival_scalps, defopt_scalp_r);


