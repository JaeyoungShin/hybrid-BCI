clear all; clc; close all;

AutobsseogDir = fullfile('C:','Program Files','MATLAB','R2013b','toolbox','eeglab12_0_2_6b','plugins','eeglab_plugin_aar-master');
MyToolboxDir = fullfile('C:','Program Files','MATLAB','R2013b','toolbox','bbci_public-master');
EegMyDataDir = fullfile('G:','IEEEdataset','rawdata','EEG');
IEEEDir = fullfile('C:','Users','Shin','Documents','MATLAB','IEEEdataset');

cd(MyToolboxDir);
startup_bbci_toolbox('DataDir',EegMyDataDir,'TmpDir','/tmp/');
cd(IEEEDir);

addpath(AutobsseogDir);

%% initial parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subdir_list = {'subject1'}
basename_list = {'imagery1','arithmetic1','imagery1','arithmetic2','imagery1','arithmetic3'};
stimDef.eeg= {16, 32;'condition1','condition2'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify discriminative frequency band. Use alpha band for example.
band_csp = [8 14]; % alpha-band

% motor imagery
for vp = 1 : length(subdir_list.imag)
    loadDir = fullfile(EegMyDataDir,subdir_list{vp});
    cd(loadDir);
    load cnt; load mrk; load mnt;
    cd(WorkingDir);
    
    cnt_temp = cnt; mrk_temp = mrk; % backup
    
    [cnt13, mrk13] = proc_appendCnt(cnt_temp{1}, cnt_temp{3}, mrk_temp{1}, mrk_temp{3});  
    eval(['[cnt_all.imag.',subdir_list.imag{vp},', mrk_all.imag.',subdir_list.imag{vp},'] = proc_appendCnt(cnt13, cnt_temp{5}, mrk13, mrk_temp{5});']);
    mnt_all = mnt;
    clear cnt mrk cnt13 mrk13
    
    [b,a]= butter(5, band_csp.imag(vp,:)/200*2);
    disp(band_csp.imag(vp,:));
    eval(['cnt_all.imag.',subdir_list.imag{vp},' = proc_filt(cnt_all.imag.', subdir_list.imag{vp},', b, a);']);
end

for vp = 1 : length(subdir_list.ment)
    disp([subdir_list.ment{vp}, ' was started']);
    loadDir = fullfile(EegMyDataDir,subdir_list.ment{vp});
    cd(loadDir);
    load cnt; load mrk; load mnt;
    cd(IEEEDir);
    cnt_temp = cnt; mrk_temp = mrk;
    
    [cnt24, mrk24] = proc_appendCnt(cnt_temp{2}, cnt_temp{4}, mrk_temp{2}, mrk_temp{4});
    eval(['[cnt_all.ment.',subdir_list.ment{vp},', mrk_all.ment.',subdir_list.ment{vp},'] = proc_appendCnt(cnt24, cnt_temp{6}, mrk24, mrk_temp{6});']);
    mnt_all = mnt;
    clear cnt mrk cnt24 mrk24
    
    [b,a]= butter(5, band_csp.ment(vp,:)/200*2);
    disp(band_csp.ment(vp,:));
    eval(['cnt_all.ment.',subdir_list.ment{vp},' = proc_filt(cnt_all.ment.', subdir_list.ment{vp},', b, a);']);
end

% initiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt_grand.imag = cnt_all.imag.VP003; mrk_grand.imag = mrk_all.imag.VP003;
cnt_grand.ment = cnt_all.ment.VP005; mrk_grand.ment = mrk_all.ment.VP005;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imagery
for vp = 2 : length(subdir_list.imag)
    disp([subdir_list.imag{vp}, ' was started']);
    eval(['[cnt_grand.imag, mrk_grand.imag] = proc_appendCnt(cnt_grand.imag, cnt_all.imag.',subdir_list.imag{vp},', mrk_grand.imag, mrk_all.imag.', subdir_list.imag{vp},');']);
end

% mental arithmetic
for vp = 2 : length(subdir_list.ment)
    disp([subdir_list.imag{vp}, ' was started']);
    eval(['[cnt_grand.ment, mrk_grand.ment] = proc_appendCnt(cnt_grand.ment, cnt_all.ment.',subdir_list.ment{vp},', mrk_grand.ment, mrk_all.ment.', subdir_list.ment{vp},');']);
end

mnt_grand = mnt_all;

ival_erd = [-8 24]*1000;
ival_scalps= [-7:5:23]*1000;


% imagery
epo= proc_segmentation(cnt_grand.imag, mrk_grand.imag, ival_erd);
epo.className = {'LMI','RMI'};
erd_car= proc_commonAverageReference(epo);
erd_car= proc_envelope(erd_car, 'MovAvgMsec', 500);
erd_car_r= proc_rSquareSigned(erd_car);

figure(4)
plot_scalpEvolution(erd_car_r, mnt_grand, ival_scalps, defopt_scalp_r,'Resolution',300, 'ContourfLevels',100);
clear epo erd_car erd_car_r

% mental arithmetic
epo= proc_segmentation(cnt_grand.ment, mrk_grand.ment, ival_erd);
epo.className = {'MA','Rest'};
erd_car= proc_commonAverageReference(epo);
erd_car= proc_envelope(erd_car, 'MovAvgMsec', 500);
erd_car_r= proc_rSquareSigned(erd_car);

figure(5)
plot_scalpEvolution(erd_car_r, mnt_grand, ival_scalps, defopt_scalp_r,'Resolution',300, 'ContourfLevels',100);



% Bandpass to the frequency band of interest
% [b,a]= butter(5, band_erd/cnt.fs*2);
% cnt= proc_filt(cnt, b, a);
%
% % Select classes 'left' and 'right'
% mrk = mrk_selectClasses(mrk,classes);
%
% % Artifact rejection based on variance criterion
% mrk= reject_varEventsAndChannels(cnt, mrk, ival_erd, ...
%                                  'DoBandpass', 0, ...
%                                  'Verbose', 1);
%
% epo= proc_segmentation(cnt, mrk, ival_erd);
% erd_lar= proc_localAverageReference(epo, mnt, 'Radius',0.4);
% erd_lar= proc_envelope(erd_lar, 'MovAvgMsec', 200);
% erd_lar= proc_baseline(erd_lar, [0 750], 'trialwise', 0);
% erd= proc_envelope(epo, 'MovAvgMsec', 200);
% erd= proc_baseline(erd, [0 750], 'trialwise', 0);
% erd_lar_r= proc_rSquareSigned(erd_lar);
% erd_r= proc_rSquareSigned(erd);
%
% fig_set(1)
% H= grid_plot(erd, mnt, defopt_erps);
% grid_addBars(erd_r, 'HScale',H.scale);
% fig_set(2)
% H= grid_plot(erd_lar, mnt, defopt_erps);
% grid_addBars(erd_lar_r, 'HScale',H.scale);
%
% fig_set(3);
% H= plot_scalpEvolutionPlusChannel(erd, mnt, clab, ival_scalps, ...
%                                   defopt_scalp_erp, ...
%                                   'ExtrapolateToMean', 1);
% grid_addBars(erd_r);
%
% fig_set(4, 'Resize',[1 2/3]);
% plot_scalpEvolution(erd_r, mnt, ival_scalps, defopt_scalp_r);
