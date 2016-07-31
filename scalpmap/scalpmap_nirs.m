% This file is for eeg-data r^2 value scalp map
% Most of MATLAB functions are available in BBCI toolbox
% Some minor code modifications might be applied
% For more tutorials, visit BBCI toolbox (https://github.com/bbci/bbci_public)

% specify your nirs data directory (NirsMyDataDir) and temporary directory (TemDir)
startup_bbci_toolbox('DataDir',EegMyDataDir,'TmpDir',TemDir);
BTB.History = 0; % to aviod error for merging cnt
fs = 10 % downsampling rate: 10 Hz

% initial parameter
subdir_list = {'subject1'}
basename_list = {'imagery1','arithmetic1','imagery2','arithmetic2','imagery3','arithmetic3'};
stimDef.eeg= {1, 2;'condition1','condition2'};

% Load data
% motor imagery
loadDir = fullfile(NirsMyDataDir,subdir_list.imag{vp});
cd(loadDir);
load cnt; load mrk; load mnt;
cd(WorkingDir);
cnt_temp = cnt; mrk_temp = mrk;

% Merge cnts
[cnt.imag.deoxy, mrk.imag.deoxy] = proc_appendCnt({cnt_temp.deoxy{1}, cnt_temp.deoxy{3}, cnt_temp.deoxy{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}});
[cnt.imag.oxy,mrk.imag.oxy]      = proc_appendCnt({cnt_temp.oxy{1}, cnt_temp.oxy{3}, cnt_temp.oxy{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}});

% Band-pass filtering
[b,a]= butter(3, [0.01 0.2]/fs*2);
cnt.imag.deoxy = proc_filt(cnt.imag.deoxy, b, a);
cnt.imag.oxy = proc_filt(cnt.imag.deoxy, b, a);

%---------------------------------------------------------------------------------------------------------------------------------------
% mental arithmetic
% Merge cnts
[cnt.ment.deoxy, mrk.ment.deoxy] = proc_appendCnt({cnt_temp.deoxy{2}, cnt_temp.deoxy{4}, cnt_temp.deoxy{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}});
[cnt.ment.oxy,mrk.ment.oxy]      = proc_appendCnt({cnt_temp.oxy{2}, cnt_temp.oxy{4}, cnt_temp.oxy{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}});

% Band-pass filtering
[b,a]= butter(3, [0.01 0.2]/fs*2);
cnt.ment.deoxy = proc_filt(cnt.ment.deoxy, b, a);
cnt.ment.oxy = proc_filt(cnt.ment.deoxy, b, a);

ival_epo = [-10 25]*1000;
ival_scalps = [-8 -5; -3 0; 2 5; 7 10; 12 15; 17 20]*1000;

% motor imagery
epo.imag.deoxy = proc_segmentation(cnt.imag.deoxy, mrk.imag.deoxy, ival_epo);
epo.imag.deoxy.className = {'LMI','RMI'};
epo_r.imag.deoxy= proc_rSquareSigned(epo.imag.deoxy);

epo.imag.oxy = proc_segmentation(cnt.imag.oxy, mrk.imag.oxy, ival_epo);
epo.imag.oxy.className = {'LMI','RMI'};
epo_r.imag.oxy= proc_rSquareSigned(epo.imag.oxy);

fig_set(1)
plot_scalpEvolution(epo.imag.deoxy, mnt, ival_scalps, defopt_scalp_r);

fig_set(2)
plot_scalpEvolution(epo.imag.oxy, mnt, ival_scalps, defopt_scalp_r);

%---------------------------------------------------------------------------------------------------------------------------------------
% mental arithmetic
epo.ment.deoxy = proc_segmentation(cnt.ment.deoxy, mrk.ment.deoxy, ival_epo);
epo.ment.deoxy.className = {'LMI','RMI'};
epo_r.ment.deoxy= proc_rSquareSigned(epo.ment.deoxy);

epo.ment.oxy = proc_segmentation(cnt.ment.oxy, mrk.ment.oxy, ival_epo);
epo.ment.oxy.className = {'LMI','RMI'};
epo_r.ment.oxy= proc_rSquareSigned(epo.ment.oxy);

fig_set(1)
plot_scalpEvolution(epo.ment.deoxy, mnt, ival_scalps, defopt_scalp_r);

fig_set(2)
plot_scalpEvolution(epo.ment.oxy, mnt, ival_scalps, defopt_scalp_r);
