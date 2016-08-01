% This file is for nirs-data classification
% Most of MATLAB functions are available in BBCI toolbox
% Some minor code modifications might be applied
% We do not guarantee all of functions works properly in your platform
% If you want to see more tutorials, visit BBCI toolbox (https://github.com/bbci/bbci_public)

% specify your eeg (EegMyDataDir) and nirs data directory (NirsMyDataDir), temporary directory (TemDir) and working directory (WorkingDir)
startup_bbci_toolbox('DataDir',NirsMyDataDir,'TmpDir','TemDir);
BTB.History = 0; % to aviod error for merging cnt

% parameters
subdir_list = {'VP001'}; % subject
basename_list = {'motor_imagery1','mental_arithmetic1','motor_imagery2','mental_arithmetic1','motor_imagery3','mental_arithmetic1'}; % task type: motor imagery / recording session: 1 - 3

stimDef.nirs = {1,2; 'condition1','condition2'};
stimDef.eeg  = {16,32; 'condition1','condition2'};

% load nirs data
loadDir = fullfile(NirsMyDataDir, subdir_list{1});
cd(loadDir);
load cnt; load mrk, load mnt; % load continous eeg signal (cnt), marker (mrk) and montage (mnt)
cd(WorkingDir)

% merge cnts in each session
% for motor imagery: imag
% for mental arithmetic: ment
cnt_temp = cnt; mrk_temp = mrk; % save data temporarily
clear cnt mrk;

[cnt.imag.deoxy, mrk.imag.deoxy] = proc_appendCnt({cnt_temp.deoxy{1}, cnt_temp.deoxy{3}, cnt_temp.deoxy{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
[cnt.imag.oxy, mrk.imag.oxy]     = proc_appendCnt({cnt_temp.oxy{1}, cnt_temp.oxy{3}, cnt_temp.oxy{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
[cnt.ment.deoxy, mrk.ment.deoxy] = proc_appendCnt({cnt_temp.deoxy{2}, cnt_temp.deoxy{4}, cnt_temp.deoxy{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}}); % merged mental arithmetic cnts
[cnt.ment.oxy, mrk.ment.oxy]     = proc_appendCnt({cnt_temp.oxy{2}, cnt_temp.oxy{4}, cnt_temp.oxy{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}}); % merged mental arithmetic cnts

% band-pass filtering
band_freq = [0.01 0.1] % Hz
ord.nirs = 3;
[z,p,k] = cheby2(ord, 30, band_freq/cnt.deoxy{session}.fs*2, 'bandpass');
[SOS,G] = zp2sos(z,p,k);

cnt.imag.deoxy = proc_filtfilt(cnt.imag.deoxy, SOS, G);
cnt.imag.oxy   = proc_filtfilt(cnt.imag.oxy, SOS, G);
cnt.ment.deoxy = proc_filtfilt(cnt.ment.deoxy, SOS, G);
cnt.ment.oxy   = proc_filtfilt(cnt.ment.oxy, SOS, G);

% channel selection
FrontalChannel = {'AF7Fp1','AF3Fp1','AF3AFz','FpzFp1','FpzAFz','FpzFp2','AF4AFz','AF4Fp2','AF8Fp2'};
MotorChannel = {'C5CP5','C5FC5','C5C3','FC3FC5','FC3C3','FC3FC1','CP3CP5','CP3C3','CP3CP1','C1C3','C1FC1','C1CP1','C2FC2','C2CP2','C2C4','FC4FC2','FC4C4','FC4FC6','CP4CP6','CP4CP2','CP4C4','C6CP6','C6C4','C6FC6'};
OccipitalChannel = {'OzPOz','OzO1','OzO2'};

cnt_org.imag.deoxy = cnt.imag.deoxy; % backup
cnt_org.imag.oxy   = cnt.imag.oxy; % backup
cnt_org.ment.deoxy = cnt.ment.deoxy; % backup
cnt_org.ment.oxy   = cnt.ment.oxy; % backup

cnt.imag.deoxy = proc_selectChannels(cnt.imag.deoxy, [MotorChannel,ParientalChannel]);
cnt.imag.oxy   = proc_selectChannels(cnt.imag.oxy, [MotorChannel,ParientalChannel]);
cnt.ment.deoxy = proc_selectChannels(cnt.ment.deoxy, [MotorChannel,ParientalChannel]);
cnt.ment.oxy   = proc_selectChannels(cnt.ment.oxy, [MotorChannel,ParientalChannel]);

% segmentation (epoching)
ival_epo  = [-10 25]*1000; % epoch range (unit: msec)
ival_base = [-3 0]*1000; % baseline correction range (unit: msec)

epo.imag.deoxy = proc_segmentation(cnt.imag.deoxy, mrk.imag.deoxy, ival_epo);
epo.imag.oxy   = proc_segmentation(cnt.imag.oxy, mrk.imag.oxy, ival_epo);
epo.imag.deoxy = proc_baseline(epo.imag.deoxy,ival_base);
epo.imag.oxy   = proc_baseline(epo.imag.oxy,ival_base);

epo.ment.deoxy = proc_segmentation(cnt.ment.deoxy, mrk.ment.deoxy, ival_epo);
epo.ment.oxy   = proc_segmentation(cnt.ment.oxy, mrk.ment.oxy, ival_epo);
epo.ment.deoxy = proc_baseline(epo.ment.deoxy,ival_base);
epo.ment.oxy   = proc_baseline(epo.ment.oxy,ival_base);

%--------------------------------------------------------------------------------------

% load occular artifact-free eeg data
loadDir = fullfile(EegMyDataDir, subdir_list{1});
cd(loadDir);
load cnt; load mrk, load mnt; % load continous eeg signal (cnt), marker (mrk) and montage (mnt)
cd(WorkingDir)

% merge cnts in each session
% for motor imagery: imag
% for mental arithmetic: ment
cnt_temp = cnt; mrk_temp = mrk; % save data temporarily
clear cnt mrk;
[cnt.imag.eeg, mrk.imag.eeg] = proc_appendCnt({cnt_temp{1}, cnt_temp{3}, cnt_temp{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
[cnt.ment.eeg, mrk.ment.eeg] = proc_appendCnt({cnt_temp{2}, cnt_temp{4}, cnt_temp{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}}); % merged mental arithmetic cnts

% Select EEG channels only (excluding EOG channels) for classification
% clab = {'F7','FAF5','F3','AFp1','AFp2','FAF6','F4','F8','FAF1','FAF2','Cz','Pz','CFC5','CFC3','CCP5','CCP3','T7','P7','P3','PPO1','OPO1','OPO2','PPO2','P4','CFC4','CFC6','CCP4','CCP6','P8','T8','VEOG','HEOG'}
cnt.imag.eeg = proc_selectChannels(cnt.imag.eeg,'not','*EOG'); % remove EOG channels (VEOG, HEOG)
mnt.imag.eeg = mnt_setElectrodePositions(cnt.imag.eeg.clab); % update montage
cnt.ment.eeg = proc_selectChannels(cnt.ment.eeg,'not','*EOG');
mnt.ment.eeg = mnt_setElectrodePositions(cnt.ment.eeg.clab);

% common average reference
cnt.imag.eeg = proc_commonAverageReference(cnt.imag.eeg);
cnt.ment.eeg = proc_commonAverageReference(cnt.ment.eeg);

% segmentation (epoching)
ival_epo  = [-10 25]*1000; % epoch range (unit: msec)
ival_base = [-3 0]*1000; % baseline correction range (unit: msec)

epo.imag.eeg = proc_segmentation(cnt.imag.eeg, mrk.imag.eeg, ival_epo);
epo.imag.eeg = proc_baseline(epo.imag.eeg,ival_base);
epo.ment.eeg = proc_segmentation(cnt.ment.eeg, mrk.ment.eeg, ival_epo);
epo.ment.eeg = proc_baseline(epo.ment.eeg,ival_base);

% frequency band selection for common spatial pattern (CSP)
MotorChannel = {'CFC5','CFC6','CFC3','CFC4','Cz,','CCP5','CCP6','CCP3','CCP4'};
ParientalChannel = {'Pz','P3','P4','P7','P8'};
FrontalChannel = {'F7','FAF5','F3','AFp1','FAF1','AFp2','FAF2','FAF6','F4','F8'};
OccipitalChannel = {'PPO1','OPO1','OPO2','PPO2'};

% channel selection
cnt_org.imag.eeg = cnt.imag.eeg; % backup
cnt.imag.eeg = proc_selectChannels(cnt.imag.eeg, [MotorChannel,ParientalChannel]);
cnt_org.ment.eeg = cnt.ment.eeg; % backup
cnt.ment.eeg = proc_selectChannels(cnt.ment.eeg, [FrontalChannel,ParientalChannel]);

% narrow frequency band selection for CSP
band_csp.imag = select_bandnarrow(cnt.imag.eeg, mrk.imag.eeg, [0 10]*1000); % band selection using 0~10 sec epoch for motor imagery
band_csp.ment = select_bandnarrow(cnt.ment.eeg, mrk.ment.eeg, [0 10]*1000); % band selection using ~10 sec epoch for mental arithmetic

% Cheby2 bandpass filter with a passband of band_csp, with at most Rp dB of passband ripple and at least Rs dB attenuation in the stopbands that are 3 Hz wide on both sides of the passband
Wp.imag = band_csp.imag/epo.imag.fs*2;
Ws.imag = [band_csp.imag(1)-3, band_csp.imag(end)+3]/epo.imag.fs*2;
Rp.imag = 3; % in dB
Rs.imag = 30; % in dB
[ord.imag, Ws.imag] = cheb2ord(Wp.imag, Ws.imag, Rp.imag, Rs.imag);
[filt_b.imag,filt_a.imag] = cheby2(ord.imag, Rs.imag, Ws.imag);
% ------------------------------------------------------------------
Wp.ment = band_csp.ment/epo.ment.fs*2;
Ws.ment = [band_csp.ment(1)-3, band_csp.ment(end)+3]/epo.ment.fs*2;
Rp.ment = 3; % in dB
Rs.ment = 30; % in dB
[ord.ment, Ws.ment] = cheb2ord(Wp.ment, Ws.ment, Rp.ment, Rs.ment);
[filt_b.ment,filt_a.ment] = cheby2(ord.ment, Rs.ment, Ws.ment);

epo.imag.eeg = proc_filtfilt(epo.imag.eeg, filt_b.imag, filt_a.imag);
epo.ment.eeg = proc_filtfilt(epo.ment.eeg, filt_b.ment, filt_a.ment);

%% classification by using moving time windows
StepSize = 1*1000; % msec
WindowSize = 3*1000; % msec
ival_start = (ival_epo(1):StepSize:ival_epo(end)-WindowSize)';
ival_end = ival_start+WindowSize;
ival = [ival_start, ival_end];
nStep = length(ival);

% average
for stepIdx = 1:nStep
    ave.imag.deoxy{stepIdx} = proc_meanAcrossTime(epo.imag.deoxy, ival(stepIdx,:));
    ave.imag.oxy{stepIdx}   = proc_meanAcrossTime(epo.imag.oxy,   ival(stepIdx,:));
    ave.ment.deoxy{stepIdx} = proc_meanAcrossTime(epo.ment.deoxy, ival(stepIdx,:));
    ave.ment.oxy{stepIdx}   = proc_meanAcrossTime(epo.ment.oxy,   ival(stepIdx,:));
end

% slope
for stepIdx = 1:nStep
    slope.imag.deoxy{stepIdx} = proc_slopeAcrossTime(epo.imag.deoxy, ival(stepIdx,:));
    slope.imag.oxy{stepIdx}   = proc_slopeAcrossTime(epo.imag.oxy,   ival(stepIdx,:));
    slope.ment.deoxy{stepIdx} = proc_slopeAcrossTime(epo.ment.deoxy, ival(stepIdx,:));
    slope.ment.oxy{stepIdx}   = proc_slopeAcrossTime(epo.ment.oxy,   ival(stepIdx,:));
end

% segment for log-variance
for stepIdx = 1:nStep
    segment.imag{stepIdx} = proc_selectIval(epo.imag.eeg, ival(stepIdx,:));
    segment.ment{stepIdx} = proc_selectIval(epo.ment.eeg, ival(stepIdx,:));
end

%--------------------------------------------------------------------------------------------------
%% 10x5-fold cross-validation
nShift = 10;
nFold = 5;

group.imag = epo.imag.deoxy.y; % epo.imag.deoxy.y == epo.imag.oxy.y == epo.imag.eeg.y
group.ment = epo.ment.deoxy.y; % epo.ment.deoxy.y == epo.ment.oxy.y

% motor imagery
for shiftIdx = 1:nShift
    indices.imag{shiftIdx} = crossvalind('Kfold',full(vec2ind(group.imag)),nFold);
    
    for stepIdx = 1:nStep
        for foldIdx = 1:nFold
            test = (indices.imag{shiftIdx} == foldIdx); train = ~test;
            
            % HbR
            x_train.deoxy.x    = [squeeze(ave.imag.deoxy{stepIdx}.x(:,:,train)); squeeze(slope.imag.deoxy{stepIdx}.x(:,:,train))];
            x_train.deoxy.y    = squeeze(ave.imag.deoxy{stepIdx}.y(:,train));
            x_train.deoxy.clab = ave.imag.deoxy{stepIdx}.clab;
            
            x_test.deoxy.x    = [squeeze(ave.imag.deoxy{stepIdx}.x(:,:,test)); squeeze(slope.imag.deoxy{stepIdx}.x(:,:,test))];
            x_test.deoxy.y    = squeeze(ave.imag.deoxy{stepIdx}.y(:,test));
            x_test.deoxy.clab = ave.imag.deoxy{stepIdx}.clab;

            % HbO
            x_train.oxy.x      = [squeeze(ave.imag.oxy{stepIdx}.x(:,:,train)); squeeze(slope.imag.oxy{stepIdx}.x(:,:,train))];
            x_train.oxy.y      = squeeze(ave.imag.oxy{stepIdx}.y(:,train));
            x_train.oxy.clab   = ave.imag.oxy{stepIdx}.clab;
         
            x_test.oxy.x    = [squeeze(ave.imag.oxy{stepIdx}.x(:,:,test)); squeeze(slope.imag.oxy{stepIdx}.x(:,:,test))];
            x_test.oxy.y    = squeeze(ave.imag.oxy{stepIdx}.y(:,test));
            x_test.oxy.clab = ave.imag.oxy{stepIdx}.clab;
            
            % eeg
            x_train.eeg.x    = squeeze(segment.imag{stepIdx}.x(:,:,train));
            x_train.eeg.y    = squeeze(segment.imag{stepIdx}.y(:,train));
            x_train.eeg.clab = segment.imag{stepIdx}.clab;
            
            x_test.eeg.x    = squeeze(segment.imag{stepIdx}.x(:,:,test));
            x_test.eeg.y    = squeeze(segment.imag{stepIdx}.y(:,test));
            x_test.eeg.clab = segment.imag{stepIdx}.clab;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% CSP for EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [csp_train, CSP_W, CSP_EIG, CSP_A] = proc_cspAuto(x_train.eeg);
            csp_train.x = csp_train.x(:,[1 2 end-1 end],:);
            
            for testIdx = 1 : size(find(test==1),1)
                csp_test.x(:,:,testIdx) = x_test.eeg.x(:,:,testIdx)*CSP_W;
            end
            csp_test.x = csp_test.x(:,[1 2 end-1 end],:);
            
            csp_train.y    = x_train.eeg.y;
            csp_train.clab = x_train.eeg.clab;
            csp_test.y     = x_test.eeg.y;
            csp_test.clab  = x_test.eeg.clab;
            
            % variance and logarithm
            var_train = proc_variance(csp_train);
            var_test  = proc_variance(csp_test);
            logvar_train = proc_logarithm(var_train);
            logvar_test  = proc_logarithm(var_test);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % feature vector
            fv_train.deoxy.x = x_train.deoxy.x; fv_train.deoxy.y = x_train.deoxy.y; fv_train.deoxy.className = {'LMI','RMI'};
            fv_test.deoxy.x  = x_test.deoxy.x;  fv_test.deoxy.y  = x_test.deoxy.y;  fv_test.deoxy.className  = {'LMI','RMI'};
            fv_train.oxy.x   = x_train.oxy.x;   fv_train.oxy.y   = x_train.oxy.y;   fv_train.oxy.className   = {'LMI','RMI'};
            fv_test.oxy.x    = x_test.oxy.x;    fv_test.oxy.y    = x_test.oxy.y;    fv_test.oxy.className    = {'LMI','RMI'};
            fv_train.eeg.x   = logvar_train.x;   fv_train.eeg.y  = x_train.eeg.y;   fv_train.eeg.className   = {'LMI','RMI'};
            fv_test.eeg.x    = logvar_test.x;   fv_test.eeg.y    = x_test.eeg.y;    fv_test.eeg.className    = {'LMI','RMI'};
            
            % for eeg only
            fv_train.eeg.x = squeeze(fv_train.eeg.x);
            fv_test.eeg.x  = squeeze(fv_test.eeg.x);
            
            y_train  = group.imag(:,train);
            y_test   = vec2ind(group.imag(:,test));
            
            % train classifier
            C.deoxy = train_RLDAshrink(fv_train.deoxy.x, y_train);
            C.oxy   = train_RLDAshrink(fv_train.oxy.x  , y_train);
            C.eeg   = train_RLDAshrink(fv_train.eeg.x  , y_train);
            
            %%%%%%%%%%%%%%%%%%%%%% train meta-classifier %%%%%%%%%%%%%%%%%%%%%%%%%
            map_train.deoxy.x = LDAmapping(C.deoxy, fv_train.deoxy.x, 'meta');
            map_train.oxy.x   = LDAmapping(C.oxy,   fv_train.oxy.x,   'meta');
            map_train.eeg.x   = LDAmapping(C.eeg,   fv_train.eeg.x,   'meta');
            
            map_test.deoxy.x  = LDAmapping(C.deoxy, fv_test.deoxy.x,  'meta');
            map_test.oxy.x    = LDAmapping(C.oxy,   fv_test.oxy.x,    'meta');
            map_test.eeg.x    = LDAmapping(C.eeg,   fv_test.eeg.x,    'meta');
            
            % meta1: HbR+HbO / meta2: HbR+EEG / meta3: HbO+EEG / meta4: HbR+HbO+EEG
            fv_train.meta1.x = [map_train.deoxy.x; map_train.oxy.x];
            fv_test.meta1.x  = [map_test.deoxy.x ; map_test.oxy.x];
            
            fv_train.meta2.x = [map_train.deoxy.x; map_train.eeg.x];
            fv_test.meta2.x  = [map_test.deoxy.x ; map_test.eeg.x];
            
            fv_train.meta3.x = [map_train.oxy.x; map_train.eeg.x];
            fv_test.meta3.x  = [map_test.oxy.x ; map_test.eeg.x];
            
            fv_train.meta4.x = [map_train.deoxy.x; map_train.oxy.x; map_train.eeg.x];
            fv_test.meta4.x  = [map_test.deoxy.x ; map_test.oxy.x ; map_test.eeg.x];
            
            y_map_train = y_train;
            y_map_test  = y_test;
            
            C.meta1 = train_RLDAshrink(fv_train.meta1.x, y_map_train);
            C.meta2 = train_RLDAshrink(fv_train.meta2.x, y_map_train);
            C.meta3 = train_RLDAshrink(fv_train.meta3.x, y_map_train);
            C.meta4 = train_RLDAshrink(fv_train.meta4.x, y_map_train);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % classification
            grouphat.deoxy(foldIdx,:)  = LDAmapping(C.deoxy,fv_test.deoxy.x);
            grouphat.oxy(foldIdx,:)    = LDAmapping(C.oxy,  fv_test.oxy.x);
            grouphat.eeg(foldIdx,:)    = LDAmapping(C.eeg,  fv_test.eeg.x);
            grouphat.meta1(foldIdx,:)  = LDAmapping(C.meta1, fv_test.meta1.x);
            grouphat.meta2(foldIdx,:)  = LDAmapping(C.meta2, fv_test.meta2.x);
            grouphat.meta3(foldIdx,:)  = LDAmapping(C.meta3, fv_test.meta3.x);
            grouphat.meta4(foldIdx,:)  = LDAmapping(C.meta4, fv_test.meta4.x);
            
            cmat.deoxy(:,:,foldIdx)  = confusionmat(y_test, grouphat.deoxy(foldIdx,:));
            cmat.oxy(:,:,foldIdx)    = confusionmat(y_test, grouphat.oxy(foldIdx,:));
            cmat.eeg(:,:,foldIdx)    = confusionmat(y_test, grouphat.eeg(foldIdx,:));
            cmat.meta1(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta1(foldIdx,:));
            cmat.meta2(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta2(foldIdx,:));
            cmat.meta3(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta3(foldIdx,:));
            cmat.meta4(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta4(foldIdx,:));
        end
        
        acc.imag.deoxy(shiftIdx,stepIdx)    = trace((sum(cmat.deoxy,3))) / sum(sum(sum(cmat.deoxy,3),2),1);
        acc.imag.oxy(shiftIdx,stepIdx)      = trace((sum(cmat.oxy,3)))   / sum(sum(sum(cmat.oxy,3),2),1);
        acc.imag.eeg(shiftIdx,stepIdx)      = trace((sum(cmat.eeg,3)))   / sum(sum(sum(cmat.eeg,3),2),1);
        acc.imag.meta1(shiftIdx,stepIdx)    = trace((sum(cmat.meta1,3)))  / sum(sum(sum(cmat.meta1,3),2),1);
        acc.imag.meta2(shiftIdx,stepIdx)    = trace((sum(cmat.meta2,3)))  / sum(sum(sum(cmat.meta2,3),2),1);
        acc.imag.meta3(shiftIdx,stepIdx)    = trace((sum(cmat.meta3,3)))  / sum(sum(sum(cmat.meta3,3),2),1);
        acc.imag.meta4(shiftIdx,stepIdx)    = trace((sum(cmat.meta4,3)))  / sum(sum(sum(cmat.meta4,3),2),1);
    end
end

mean_acc.imag.deoxy  = mean(acc.imag.deoxy,1)';
mean_acc.imag.oxy    = mean(acc.imag.oxy,1)';
mean_acc.imag.eeg    = mean(acc.imag.eeg,1)';
mean_acc.imag.meta1  = mean(acc.imag.meta1,1)';
mean_acc.imag.meta2  = mean(acc.imag.meta2,1)';
mean_acc.imag.meta3  = mean(acc.imag.meta3,1)';
mean_acc.imag.meta4  = mean(acc.imag.meta4,1)';
    
time = (ival(:,2)/1000)';
acc.imag.time = time;

% display
figure('Name','Hybrid Motor imagery', 'Number', 'off')
subplot(2,2,1); plot(time, mean_acc.imag.deoxy, 'r', time, mean_acc.imag.oxy, 'b', time, mean_acc.imag.eeg,'k','Location','NorthWest');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MI DEOXY','MI OXY','MI EEG','Location','NorthWest');
subplot(2,2,2); plot(time, mean_acc.imag.deoxy, 'r', time, mean_acc.imag.meta2, 'b');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MI DEOXY','MI DEOXY+EEG','Location','NorthWest');
subplot(2,2,3); plot(time, mean_acc.imag.oxy,   'r', time, mean_acc.imag.meta3, 'b');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MI OXY','MI OXY+EEG','Location','NorthWest');
subplot(2,2,4); plot(time, mean_acc.imag.meta1, 'r', time, mean_acc.imag.meta4, 'b');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MI DEOXY+OXY','MI DEOXY+OXY+EEG','Location','NorthWest');

%-----------------------------------------------------------------------------------------------------------
% mental arithmetic
for shiftIdx = 1:nShift
    indices.ment{shiftIdx} = crossvalind('Kfold',full(vec2ind(group.ment)),nFold);
    
    for stepIdx = 1:nStep
        for foldIdx = 1:nFold
            test = (indices.ment{shiftIdx} == foldIdx); train = ~test;
            
            % HbR
            x_train.deoxy.x    = [squeeze(ave.ment.deoxy{stepIdx}.x(:,:,train)); squeeze(slope.ment.deoxy{stepIdx}.x(:,:,train))];
            x_train.deoxy.y    = squeeze(ave.ment.deoxy{stepIdx}.y(:,train));
            x_train.deoxy.clab = ave.ment.deoxy{stepIdx}.clab;
            
            x_test.deoxy.x    = [squeeze(ave.ment.deoxy{stepIdx}.x(:,:,test)); squeeze(slope.ment.deoxy{stepIdx}.x(:,:,test))];
            x_test.deoxy.y    = squeeze(ave.ment.deoxy{stepIdx}.y(:,test));
            x_test.deoxy.clab = ave.ment.deoxy{stepIdx}.clab;
            
            % HbO
            x_train.oxy.x      = [squeeze(ave.ment.oxy{stepIdx}.x(:,:,train)); squeeze(slope.ment.oxy{stepIdx}.x(:,:,train))];
            x_train.oxy.y      = squeeze(ave.ment.oxy{stepIdx}.y(:,train));
            x_train.oxy.clab   = ave.ment.oxy{stepIdx}.clab;
            
            x_test.oxy.x    = [squeeze(ave.ment.oxy{stepIdx}.x(:,:,test)); squeeze(slope.ment.oxy{stepIdx}.x(:,:,test))];
            x_test.oxy.y    = squeeze(ave.ment.oxy{stepIdx}.y(:,test));
            x_test.oxy.clab = ave.ment.oxy{stepIdx}.clab;
            
            % eeg
            x_train.eeg.x    = squeeze(segment.ment{stepIdx}.x(:,:,train));
            x_train.eeg.y    = squeeze(segment.ment{stepIdx}.y(:,train));
            x_train.eeg.clab = segment.ment{stepIdx}.clab;
            
            x_test.eeg.x    = squeeze(segment.ment{stepIdx}.x(:,:,test));
            x_test.eeg.y    = squeeze(segment.ment{stepIdx}.y(:,test));
            x_test.eeg.clab = segment.ment{stepIdx}.clab;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% CSP for EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [csp_train, CSP_W, CSP_EIG, CSP_A] = proc_cspAuto(x_train.eeg);
            csp_train.x = csp_train.x(:,[1 2 end-1 end],:);
            
            for testIdx = 1 : size(find(test==1),1)
                csp_test.x(:,:,testIdx) = x_test.eeg.x(:,:,testIdx)*CSP_W;
            end
            csp_test.x = csp_test.x(:,[1 2 end-1 end],:);
            
            csp_train.y    = x_train.eeg.y;
            csp_train.clab = x_train.eeg.clab;
            csp_test.y     = x_test.eeg.y;
            csp_test.clab  = x_test.eeg.clab;
            
            % variance and logarithm
            var_train = proc_variance(csp_train);
            var_test  = proc_variance(csp_test);
            logvar_train = proc_logarithm(var_train);
            logvar_test  = proc_logarithm(var_test);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % feature vector
            fv_train.deoxy.x = x_train.deoxy.x; fv_train.deoxy.y = x_train.deoxy.y; fv_train.deoxy.className = {'LMI','RMI'};
            fv_test.deoxy.x  = x_test.deoxy.x;  fv_test.deoxy.y  = x_test.deoxy.y;  fv_test.deoxy.className  = {'LMI','RMI'};
            fv_train.oxy.x   = x_train.oxy.x;   fv_train.oxy.y   = x_train.oxy.y;   fv_train.oxy.className   = {'LMI','RMI'};
            fv_test.oxy.x    = x_test.oxy.x;    fv_test.oxy.y    = x_test.oxy.y;    fv_test.oxy.className    = {'LMI','RMI'};
            fv_train.eeg.x   = logvar_train.x;   fv_train.eeg.y  = x_train.eeg.y;   fv_train.eeg.className   = {'LMI','RMI'};
            fv_test.eeg.x    = logvar_test.x;   fv_test.eeg.y    = x_test.eeg.y;    fv_test.eeg.className    = {'LMI','RMI'};
            
            % for eeg only
            fv_train.eeg.x = squeeze(fv_train.eeg.x);
            fv_test.eeg.x  = squeeze(fv_test.eeg.x);
            
            y_train  = group.ment(:,train);
            y_test   = vec2ind(group.ment(:,test));
            
            % train classifier
            C.deoxy = train_RLDAshrink(fv_train.deoxy.x, y_train);
            C.oxy   = train_RLDAshrink(fv_train.oxy.x  , y_train);
            C.eeg   = train_RLDAshrink(fv_train.eeg.x  , y_train);
            
            %%%%%%%%%%%%%%%%%%%%%% train meta-classifier %%%%%%%%%%%%%%%%%%%%%%%%%
            map_train.deoxy.x = LDAmapping(C.deoxy, fv_train.deoxy.x, 'meta');
            map_train.oxy.x   = LDAmapping(C.oxy,   fv_train.oxy.x,   'meta');
            map_train.eeg.x   = LDAmapping(C.eeg,   fv_train.eeg.x,   'meta');
            
            map_test.deoxy.x  = LDAmapping(C.deoxy, fv_test.deoxy.x,  'meta');
            map_test.oxy.x    = LDAmapping(C.oxy,   fv_test.oxy.x,    'meta');
            map_test.eeg.x    = LDAmapping(C.eeg,   fv_test.eeg.x,    'meta');
            
            % meta1: HbR+HbO / meta2: HbR+EEG / meta3: HbO+EEG / meta4: HbR+HbO+EEG
            fv_train.meta1.x = [map_train.deoxy.x; map_train.oxy.x];
            fv_test.meta1.x  = [map_test.deoxy.x ; map_test.oxy.x];
            
            fv_train.meta2.x = [map_train.deoxy.x; map_train.eeg.x];
            fv_test.meta2.x  = [map_test.deoxy.x ; map_test.eeg.x];
            
            fv_train.meta3.x = [map_train.oxy.x; map_train.eeg.x];
            fv_test.meta3.x  = [map_test.oxy.x ; map_test.eeg.x];
            
            fv_train.meta4.x = [map_train.deoxy.x; map_train.oxy.x; map_train.eeg.x];
            fv_test.meta4.x  = [map_test.deoxy.x ; map_test.oxy.x ; map_test.eeg.x];
            
            y_map_train = y_train;
            y_map_test  = y_test;
            
            C.meta1 = train_RLDAshrink(fv_train.meta1.x, y_map_train);
            C.meta2 = train_RLDAshrink(fv_train.meta2.x, y_map_train);
            C.meta3 = train_RLDAshrink(fv_train.meta3.x, y_map_train);
            C.meta4 = train_RLDAshrink(fv_train.meta4.x, y_map_train);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % classification
            grouphat.deoxy(foldIdx,:)  = LDAmapping(C.deoxy,fv_test.deoxy.x);
            grouphat.oxy(foldIdx,:)    = LDAmapping(C.oxy,  fv_test.oxy.x);
            grouphat.eeg(foldIdx,:)    = LDAmapping(C.eeg,  fv_test.eeg.x);
            grouphat.meta1(foldIdx,:)  = LDAmapping(C.meta1, fv_test.meta1.x);
            grouphat.meta2(foldIdx,:)  = LDAmapping(C.meta2, fv_test.meta2.x);
            grouphat.meta3(foldIdx,:)  = LDAmapping(C.meta3, fv_test.meta3.x);
            grouphat.meta4(foldIdx,:)  = LDAmapping(C.meta4, fv_test.meta4.x);
            
            cmat.deoxy(:,:,foldIdx)  = confusionmat(y_test, grouphat.deoxy(foldIdx,:));
            cmat.oxy(:,:,foldIdx)    = confusionmat(y_test, grouphat.oxy(foldIdx,:));
            cmat.eeg(:,:,foldIdx)    = confusionmat(y_test, grouphat.eeg(foldIdx,:));
            cmat.meta1(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta1(foldIdx,:));
            cmat.meta2(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta2(foldIdx,:));
            cmat.meta3(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta3(foldIdx,:));
            cmat.meta4(:,:,foldIdx)  = confusionmat(y_test, grouphat.meta4(foldIdx,:));
        end
        
        acc.ment.deoxy(shiftIdx,stepIdx)    = trace((sum(cmat.deoxy,3))) / sum(sum(sum(cmat.deoxy,3),2),1);
        acc.ment.oxy(shiftIdx,stepIdx)      = trace((sum(cmat.oxy,3)))   / sum(sum(sum(cmat.oxy,3),2),1);
        acc.ment.eeg(shiftIdx,stepIdx)      = trace((sum(cmat.eeg,3)))   / sum(sum(sum(cmat.eeg,3),2),1);
        acc.ment.meta1(shiftIdx,stepIdx)    = trace((sum(cmat.meta1,3)))  / sum(sum(sum(cmat.meta1,3),2),1);
        acc.ment.meta2(shiftIdx,stepIdx)    = trace((sum(cmat.meta2,3)))  / sum(sum(sum(cmat.meta2,3),2),1);
        acc.ment.meta3(shiftIdx,stepIdx)    = trace((sum(cmat.meta3,3)))  / sum(sum(sum(cmat.meta3,3),2),1);
        acc.ment.meta4(shiftIdx,stepIdx)    = trace((sum(cmat.meta4,3)))  / sum(sum(sum(cmat.meta4,3),2),1);
    end
end

mean_acc.ment.deoxy  = mean(acc.ment.deoxy,1)';
mean_acc.ment.oxy    = mean(acc.ment.oxy,1)';
mean_acc.ment.eeg    = mean(acc.ment.eeg,1)';
mean_acc.ment.meta1  = mean(acc.ment.meta1,1)';
mean_acc.ment.meta2  = mean(acc.ment.meta2,1)';
mean_acc.ment.meta3  = mean(acc.ment.meta3,1)';
mean_acc.ment.meta4  = mean(acc.ment.meta4,1)';

time = (ival(:,2)/1000)';
acc.ment.time = time;

% display
figure('Name','Hybrid Mental arithmetic', 'Number', 'off')
subplot(2,2,1); plot(time, mean_acc.ment.deoxy, 'r', time, mean_acc.ment.oxy, 'b', time, mean_acc.ment.eeg,'k','Location','NorthWest');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MA DEOXY','MA OXY','MA EEG','Location','NorthWest');
subplot(2,2,2); plot(time, mean_acc.ment.deoxy, 'r', time, mean_acc.ment.meta2, 'b');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MA DEOXY','MA DEOXY+EEG','Location','NorthWest');
subplot(2,2,3); plot(time, mean_acc.ment.oxy,   'r', time, mean_acc.ment.meta3, 'b');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MA OXY','MA OXY+EEG','Location','NorthWest');
subplot(2,2,4); plot(time, mean_acc.ment.meta1, 'r', time, mean_acc.ment.meta4, 'b');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on; legend('MA DEOXY+OXY','MA DEOXY+OXY+EEG','Location','NorthWest');
