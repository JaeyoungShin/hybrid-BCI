% This file is for applying the modified Beer-Lambert law (MBLL) to NIRS data
% Most of MATLAB functions are available in BBCI toolbox
% Some minor code modifications might be applied
% We do not guarantee all of functions works properly in your platform
% If you want to see more tutorials, visit BBCI toolbox (https://github.com/bbci/bbci_public)
% Specify your nirs data directory (NirsMyDataDir), temporary directory (TemDir) and working directory (WorkingDir)
%
% For example:
% NirsMyDataDir = fullfile('H:','DEMOdataset','rawdata','NIRS');
% TemDir = fullfile('H:','DEMOdataset','temporary');
% WorkingDir = fullfile('C:','Users','shin','Documents','MATLAB','IEEE_tutorial');
%
% 1. load cnt mrk files
% 2. Apply MBLL and save post-MBLL data
% 3. Load data and divide data into deoxy- and oxy-hemoglobin data
% 4. Merge cnts for motor imagery and mental arithmetic separately.

% load toolbox
startup_bbci_toolbox('DataDir',NirsMyDataDir,'TmpDir',TemDir, 'RawDir', NirsMyDataDir, 'MatDir', NirsMyDataDir, 'History', 0); % History must be 0 to avoid a potential error in merging cnts.

% parameters
subdir_list = {'subject 01'}; % subject
basename_list = {'motor_imagery1','mental_arithmetic1','motor_imagery2','mental_arithmetic2','motor_imagery3','mental_arithmetic3'};
stimDef.eeg = {1,2; 'condition1','condition2'};

% load NIRS data
loadDir = fullfile(NirsMyDataDir, subdir_list{1});
cd(loadDir);
load cnt; load mrk; % load continous eeg signal (cnt), marker (mrk) and montage (mnt)
cd(WorkingDir)

% MBLL and save post-MBLL data
for idx = 1 : 6
    filename{idx} = fullfile(subdir_list{1}, ['session',num2str(idx)]);
    cntHb{idx} = proc_BeerLambert(cnt{idx}, 'Opdist', 30, 'DPF', [5.98 7.15], 'Citation', 1);
    file_saveNIRSMatlab(filename{idx}, cntHb{idx}, mrk{idx}, mnt);
end

clear cnt cntHb mrk;

% load deoxy- and oxy-hemoglobin data
for idx = 1 : 6
    [cnt_temp.deoxy{idx}, mrk_temp{idx}, mnt] = file_loadNIRSMatlab(filename{idx}, 'Signal','deoxy');
    [cnt_temp.oxy{idx}  , ~, ~]   = file_loadNIRSMatlab(filename{idx}, 'Signal','oxy');
end

% merge cnts for motor imagery and mental arithmetic separately
[cnt.imag.deoxy, mrk.imag] = proc_appendCnt({cnt_temp.deoxy{1}, cnt_temp.deoxy{3}, cnt_temp.deoxy{5}},{mrk_temp{1},mrk_temp{3},mrk_temp{5}});
[cnt.imag.oxy,   mrk.imag] = proc_appendCnt({cnt_temp.oxy{1},   cnt_temp.oxy{3},   cnt_temp.oxy{5}},  {mrk_temp{1},mrk_temp{3},mrk_temp{5}});
[cnt.ment.deoxy, mrk.ment] = proc_appendCnt({cnt_temp.deoxy{2}, cnt_temp.deoxy{4}, cnt_temp.deoxy{6}},{mrk_temp{2},mrk_temp{4},mrk_temp{6}});
[cnt.ment.oxy,   mrk.ment] = proc_appendCnt({cnt_temp.oxy{2},   cnt_temp.oxy{4},   cnt_temp.oxy{6}},  {mrk_temp{2},mrk_temp{4},mrk_temp{6}});

clear cnt_temp
