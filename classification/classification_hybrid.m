% This file is for nirs-data classification
% Most of MATLAB functions are available in BBCI toolbox
% Some minor code modifications might be applied
% We do not guarantee all of functions works properly in your platform
% If you want to see more tutorials, visit BBCI toolbox (https://github.com/bbci/bbci_public)

% specify your nirs data directory (NirsMyDataDir), temporary directory (TemDir) and working directory (WorkingDir)
startup_bbci_toolbox('DataDir',NirsMyDataDir,'TmpDir','TemDir);
BTB.History = 0; % to aviod error for merging cnt
% parameters
subdir_list = {'VP001'}; % subject
basename_list = {'motor_imagery1','mental_arithmetic1','motor_imagery2','mental_arithmetic1','motor_imagery3','mental_arithmetic1'}; % task type: motor imagery / recording session: 1 - 3
stimDef.eeg = {1,2; 'condition1','condition2'};

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
ord = 3;
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

% cross-validation (nShift x nFold-fold cross-validation)
% for more convenient use, please refer to 'https://github.com/bbci/bbci_public/blob/master/demos/demo_validation_csp.m'
% cross-validation below was written for meta-classification for EEG-NIRS hybrid BCI
group.imag = epo.imag.deoxy.y; % epo.imag.deoxy.y == epo.imag.oxy.y
group.ment = epo.ment.deoxy.y; % epo.ment.deoxy.y == epo.ment.oxy.y
    
% motor imagery
for shiftIdx = 1:nShift
    indices.imag{shiftIdx} = crossvalind('Kfold',full(vec2ind(group.imag)),nFold);
    for stepIdx = 1:nStep
        fprintf('Motor imagery, Repeat: %d/%d, Step: %d/%d\n',shiftIdx, nShift, stepIdx, nStep);
        for foldIdx = 1:nFold
            test = (indices.imag{shiftIdx} == foldIdx); train = ~test;
            
            x_train.deoxy.x = [squeeze(ave.imag.deoxy{stepIdx}.x(:,:,train)); squeeze(slope.imag.deoxy{stepIdx}.x(:,:,train))];
            x_train.deoxy.y = squeeze(ave.imag.deoxy{stepIdx}.y(:,train));
            x_train.deoxy.clab = ave.imag.deoxy{stepIdx}.clab;
            
            x_train.oxy.x = [squeeze(ave.imag.oxy{stepIdx}.x(:,:,train)); squeeze(slope.imag.oxy{stepIdx}.x(:,:,train))];
            x_train.oxy.y   = squeeze(ave.imag.oxy{stepIdx}.y(:,train));
            x_train.oxy.clab   = ave.imag.oxy{stepIdx}.clab;
            
            x_test.deoxy.x = [squeeze(ave.imag.deoxy{stepIdx}.x(:,:,test)); squeeze(slope.imag.deoxy{stepIdx}.x(:,:,test))];
            x_test.deoxy.y = squeeze(ave.imag.deoxy{stepIdx}.y(:,test));
            x_test.deoxy.clab = ave.imag.deoxy{stepIdx}.clab;
            
            x_test.oxy.x = [squeeze(ave.imag.oxy{stepIdx}.x(:,:,test)); squeeze(slope.imag.oxy{stepIdx}.x(:,:,test))];
            x_test.oxy.y = squeeze(ave.imag.oxy{stepIdx}.y(:,test));
            x_test.oxy.clab = ave.imag.oxy{stepIdx}.clab;
            
            % feature vector
            fv_train.deoxy.x = x_train.deoxy.x; fv_train.deoxy.y = x_train.deoxy.y; fv_train.deoxy.className = {'Cond1','Cond2'};
            fv_test.deoxy.x  = x_test.deoxy.x;  fv_test.deoxy.y  = x_test.deoxy.y;  fv_test.deoxy.className  = {'Cond1','Cond2'};
            fv_train.oxy.x   = x_train.oxy.x;   fv_train.oxy.y   = x_train.oxy.y;   fv_train.oxy.className   = {'Cond1','Cond2'};
            fv_test.oxy.x    = x_test.oxy.x;    fv_test.oxy.y    = x_test.oxy.y;    fv_test.oxy.className    = {'Cond1','Cond2'};
            
            y_train  = group.imag(:,train);
            y_test   = vec2ind(group.imag(:,test));
            
            
            % train classifier
            C.deoxy = train_RLDAshrink(fv_train.deoxy.x,y_train);
            C.oxy   = train_RLDAshrink(fv_train.oxy.x  ,y_train);
            
            % classification
            grouphat.deoxy(foldIdx,:) = LDAmapping(C.deoxy,fv_test.deoxy.x);
            grouphat.oxy(foldIdx,:)   = LDAmapping(C.oxy,  fv_test.oxy.x);
                
            cmat.deoxy(:,:,foldIdx) = confusionmat(y_test, grouphat.deoxy(foldIdx,:));
            cmat.oxy(:,:,foldIdx)   = confusionmat(y_test, grouphat.oxy(foldIdx,:));
            
        end
        acc.deoxy.imag(shiftIdx,stepIdx)   = trace((sum(cmat.deoxy,3))) / sum(sum(sum(cmat.deoxy,3),2),1);
        acc.oxy.imag(shiftIdx,stepIdx)     = trace((sum(cmat.oxy,3))) / sum(sum(sum(cmat.oxy,3),2),1);
    end
end

mean_acc.imag.deoxy = mean(acc.deoxy.imag,1)';
mean_acc.imag.oxy = mean(acc.oxy.imag,1)';

%-----------------------------------------------------------------------------------------------------------
% mental arithmetic
for shiftIdx = 1:nShift
    indices.ment{shiftIdx} = crossvalind('Kfold',full(vec2ind(group.ment)),nFold);
    for stepIdx = 1:nStep
        fprintf('Mental arithmetic, Repeat: %d/%d, Step: %d/%d\n',shiftIdx, nShift, stepIdx, nStep);
        for foldIdx = 1:nFold
            test = (indices.ment{shiftIdx} == foldIdx); train = ~test;
            
            x_train.deoxy.x = [squeeze(ave.ment.deoxy{stepIdx}.x(:,:,train)); squeeze(slope.ment.deoxy{stepIdx}.x(:,:,train))];
            x_train.deoxy.y = squeeze(ave.ment.deoxy{stepIdx}.y(:,train));
            x_train.deoxy.clab = ave.ment.deoxy{stepIdx}.clab;
            
            x_train.oxy.x = [squeeze(ave.ment.oxy{stepIdx}.x(:,:,train)); squeeze(slope.ment.oxy{stepIdx}.x(:,:,train))];
            x_train.oxy.y   = squeeze(ave.ment.oxy{stepIdx}.y(:,train));
            x_train.oxy.clab   = ave.ment.oxy{stepIdx}.clab;
            
            x_test.deoxy.x = [squeeze(ave.ment.deoxy{stepIdx}.x(:,:,test)); squeeze(slope.ment.deoxy{stepIdx}.x(:,:,test))];
            x_test.deoxy.y = squeeze(ave.ment.deoxy{stepIdx}.y(:,test));
            x_test.deoxy.clab = ave.ment.deoxy{stepIdx}.clab;
            
            x_test.oxy.x = [squeeze(ave.ment.oxy{stepIdx}.x(:,:,test)); squeeze(slope.ment.oxy{stepIdx}.x(:,:,test))];
            x_test.oxy.y = squeeze(ave.ment.oxy{stepIdx}.y(:,test));
            x_test.oxy.clab = ave.ment.oxy{stepIdx}.clab;

            % feature vector
            fv_train.deoxy.x = x_train.deoxy.x; fv_train.deoxy.y = x_train.deoxy.y; fv_train.deoxy.className = {'Cond1','Cond2'};
            fv_test.deoxy.x  = x_test.deoxy.x;  fv_test.deoxy.y  = x_test.deoxy.y;  fv_test.deoxy.className  = {'Cond1','Cond2'};
            fv_train.oxy.x   = x_train.oxy.x;   fv_train.oxy.y   = x_train.oxy.y;   fv_train.oxy.className   = {'Cond1','Cond2'};
            fv_test.oxy.x    = x_test.oxy.x;    fv_test.oxy.y    = x_test.oxy.y;    fv_test.oxy.className    = {'Cond1','Cond2'};
            
            y_train  = group.ment(:,train);
            y_test   = vec2ind(group.ment(:,test));
            
            % train classifier
            C.deoxy = train_RLDAshrink(fv_train.deoxy.x,y_train);
            C.oxy   = train_RLDAshrink(fv_train.oxy.x  ,y_train);
            
            % classification
            grouphat.deoxy(foldIdx,:) = LDAmapping(C.deoxy,fv_test.deoxy.x);
            grouphat.oxy(foldIdx,:)   = LDAmapping(C.oxy,  fv_test.oxy.x);
            
            cmat.deoxy(:,:,foldIdx) = confusionmat(y_test, grouphat.deoxy(foldIdx,:));
            cmat.oxy(:,:,foldIdx)   = confusionmat(y_test, grouphat.oxy(foldIdx,:));
            
        end
        acc.deoxy.ment(shiftIdx,stepIdx)   = trace((sum(cmat.deoxy,3))) / sum(sum(sum(cmat.deoxy,3),2),1);
        acc.oxy.ment(shiftIdx,stepIdx)     = trace((sum(cmat.oxy,3))) / sum(sum(sum(cmat.oxy,3),2),1);
    end
end

mean_acc.ment.deoxy = mean(acc.deoxy.ment,1)';
mean_acc.ment.oxy = mean(acc.oxy.ment,1)';

result= [mean_acc.imag.deoxy, mean_acc.imag.oxy, mean_acc.ment.deoxy, mean_acc.ment.oxy];

% display
time = (ival(:,2)/1000)';
figure('Name',subdir_list{vp},'Number','off')
plot(time,result);
legend('imag deoxy','imag oxy','ment deoxy','ment oxy');
xlim([time(1) time(end)]); ylim([0.4 1]); grid on;

