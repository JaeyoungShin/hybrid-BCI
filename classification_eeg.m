% This file is for eeg-data classification
% All MATLAB functions are available in BBCI toolbox
% Some minor code modifications might be applied
% We do not guarantee all of functions works properly in your platform
% If you want to see more tutorials, visit BBCI toolbox (https://github.com/bbci/bbci_public)

% specify your eeg data directory (EegMyDataDir) and temporary directory (TemDir)
startup_bbci_toolbox('DataDir',EegMyDataDir,'TmpDir','/tmp/');
BTB.History = 0; % to aviod error for merging cnt

% parameters
subdir_list = {'VP001'}; % subject
basename_list = {'motor_imagery1','mental_arithmetic1','motor_imagery2','mental_arithmetic1','motor_imagery3','mental_arithmetic1'}; % task type: motor imagery / recording session: 1 - 3
stimDef.eeg = {16,32; 'condition1','condition2'};

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
[cnt.imag, mrk.imag] = proc_appendCnt({cnt_temp{1}, cnt_temp{3}, cnt_temp{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
[cnt.ment, mrk.ment] = proc_appendCnt({cnt_temp{2}, cnt_temp{4}, cnt_temp{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}}); % merged mental arithmetic cnts

% Select EEG channels only (excluding EOG channels) for classification
% clab = {'F7','FAF5','F3','AFp1','AFp2','FAF6','F4','F8','FAF1','FAF2','Cz','Pz','CFC5','CFC3','CCP5','CCP3','T7','P7','P3','PPO1','OPO1','OPO2','PPO2','P4','CFC4','CFC6','CCP4','CCP6','P8','T8','VEOG','HEOG'}
cnt.imag = proc_selectChannels(cnt.imag,'not','*EOG'); % remove EOG channels (VEOG, HEOG)
mnt.imag = mnt_setElectrodePositions(cnt.imag.clab); % update montage
cnt.ment = proc_selectChannels(cnt.ment,'not','*EOG');
mnt.ment = mnt_setElectrodePositions(cnt.ment.clab);

% common average reference
cnt.imag = proc_commonAverageReference(cnt.imag);
cnt.ment = proc_commonAverageReference(cnt.ment);

% segmentation (epoching)
ival_epo  = [-10 25]*1000; % epoch range (unit: msec)
ival_base = [-3 0]*1000; % baseline correction range (unit: msec)

epo.imag = proc_segmentation(cnt.imag, mrk.imag, ival_epo);
epo.imag = proc_baseline(epo.imag,ival_base);
epo.ment = proc_segmentation(cnt.ment, mrk.ment, ival_epo);
epo.ment = proc_baseline(epo.ment,ival_base);

% frequency band selection for common spatial pattern (CSP)
MotorChannel = {'CFC5','CFC6','CFC3','CFC4','Cz,','CCP5','CCP6','CCP3','CCP4'};
ParientalChannel = {'Pz','P3','P4','P7','P8'};
FrontalChannel = {'F7','FAF5','F3','AFp1','FAF1','AFp2','FAF2','FAF6','F4','F8'};
OccipitalChannel = {'PPO1','OPO1','OPO2','PPO2'};

% channel selection
cnt_org.imag = cnt.imag; % backup
cnt.imag = proc_selectChannels(cnt.imag, [MotorChannel,ParientalChannel]);
cnt_org.ment = cnt.ment; % backup
cnt.ment = proc_selectChannels(cnt.ment, [FrontalChannel,ParientalChannel]);

% narrow frequency band selection for CSP
band_csp.imag = select_bandnarrow(cnt.imag, mrk.imag, [0 10]*1000); % band selection using 0~10 sec epoch for motor imagery
band_csp.ment = select_bandnarrow(cnt.ment, mrk.ment, [0 10]*1000); % band selection using ~10 sec epoch for mental arithmetic

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

epo.imag = proc_filtfilt(epo.imag, filt_b.imag, filt_a.imag);
epo.ment = proc_filtfilt(epo.ment, filt_b.ment, filt_a.ment);

%% classification by using moving time windows
StepSize = 1*1000; % msec
WindowSize = 3*1000; % msec
ival_start = (ival_epo(1):StepSize:ival_epo(end)-WindowSize)';
ival_end = ival_start+WindowSize;
ival = [ival_start, ival_end];
nStep = length(ival);
  
% select time interval for moving time window
for stepIdx = 1:nStep
    segment.imag{stepIdx} = proc_selectIval(epo.imag, ival(stepIdx,:));
    segment.ment{stepIdx} = proc_selectIval(epo.ment, ival(stepIdx,:));
end

% cross-validation (nShift x nFold-fold cross-validation)
% for more convenient use, please refer to 'https://github.com/bbci/bbci_public/blob/master/demos/demo_validation_csp.m'
% cross-validation below was written for meta-classification for EEG-NIRS hybrid BCI
nShift = 5; % number of repeatitions
nFold = 5;  % number of subsets
group.imag = epo.imag.y;
group.ment = epo.ment.y;

% motor imagery
for shiftIdx = 1:nShift
    indices.imag{shiftIdx} = crossvalind('Kfold',full(vec2ind(group.imag)),nFold);
    for stepIdx = 1:nStep
        fprintf('Motor imagery, Repeat: %d/%d, Step: %d/%d\n',shiftIdx, nShift, stepIdx, nStep);
        for foldIdx = 1:nFold
            clear csp_train csp_test
            test = (indices.imag{shiftIdx} == foldIdx); train = ~test;
            
            x_train.x = segment.imag{stepIdx}.x(:,:,train);
            x_train.y = segment.imag{stepIdx}.y(:,train);
            x_train.clab = segment.imag{stepIdx}.clab;
            
            x_test.x = segment.imag{stepIdx}.x(:,:,test);
            x_test.y = segment.imag{stepIdx}.y(:,test);
            x_test.clab = segment.imag{stepIdx}.clab;
            
            % CSP
            [csp_train, CSP_W, CSP_EIG, CSP_A] = proc_cspAuto(x_train);
            csp_train.x = csp_train.x(:,[1 2 end-1 end],:); % select the first 2 and the last 2 CSP components
            
            for testIdx = 1 : size(find(test==1),1)
                csp_test.x(:,:,testIdx) = x_test.x(:,:,testIdx)*CSP_W;
            end
            csp_test.x = csp_test.x(:,[1 2 end-1 end],:); % select the first 2 and the last 2 CSP components
            
            csp_train.y    = x_train.y;
            csp_train.clab = x_train.clab;
            csp_test.y     = x_test.y;
            csp_test.clab  = x_test.clab;
            
            % variance and logarithm
            var_train = proc_variance(csp_train);
            var_test  = proc_variance(csp_test);
            logvar_train = proc_logarithm(var_train);
            logvar_test  = proc_logarithm(var_test);
            
            % feature vector
            fv_train.eeg = squeeze(logvar_train.x);
            fv_test.eeg = squeeze(logvar_test.x);
            
            y_train  = group.imag(:,train);
            y_test   = vec2ind(group.imag(:,test));
            
            % train classifier
            C.eeg = train_RLDAshrink(fv_train.eeg,y_train); % shrinkage linear discriminant analysis
            
            % classification
            grouphat.eeg(foldIdx,:) = LDAmapping(C.eeg,fv_test.eeg); % use custom function
            cmat.eeg(:,:,foldIdx) = confusionmat(y_test, grouphat.eeg(foldIdx,:));
        end
        
        acc.eeg.imag(shiftIdx,stepIdx)   = trace((sum(cmat.eeg,3))) / sum(sum(sum(cmat.eeg,3),2),1);
    end
end

mean_acc.eeg.imag = mean(acc.eeg.imag,1)';
% -----------------------------------------------------------------------------------------------------
% mental arithmetic
for shiftIdx = 1:nShift
    indices.ment{shiftIdx} = crossvalind('Kfold',full(vec2ind(group.ment)),nFold);
    for stepIdx = 1:nStep
        fprintf('Mental arithmetic, Repeat: %d/%d, Step: %d/%d\n',shiftIdx, nShift, stepIdx, nStep);
        for foldIdx = 1:nFold
            clear csp_train csp_test
            test = (indices.ment{shiftIdx} == foldIdx); train = ~test;
            
            x_train.x = segment.ment{stepIdx}.x(:,:,train);
            x_train.y = segment.ment{stepIdx}.y(:,train);
            x_train.clab = segment.ment{stepIdx}.clab;
            
            x_test.x = segment.ment{stepIdx}.x(:,:,test);
            x_test.y = segment.ment{stepIdx}.y(:,test);
            x_test.clab = segment.ment{stepIdx}.clab;
            
            % CSP
            [csp_train, CSP_W, CSP_EIG, CSP_A] = proc_cspAuto(x_train);
            csp_train.x = csp_train.x(:,[1 2 end-1 end],:); % select the first 2 and the last 2 CSP components
            
            for testIdx = 1 : size(find(test==1),1)
                csp_test.x(:,:,testIdx) = x_test.x(:,:,testIdx)*CSP_W;
            end
            csp_test.x = csp_test.x(:,[1 2 end-1 end],:); % select the first 2 and the last 2 CSP components
            
            csp_train.y    = x_train.y;
            csp_train.clab = x_train.clab;
            csp_test.y     = x_test.y;
            csp_test.clab  = x_test.clab;
            
            % variance and logarithm
            var_train = proc_variance(csp_train);
            var_test  = proc_variance(csp_test);
            logvar_train = proc_logarithm(var_train);
            logvar_test  = proc_logarithm(var_test);
            
            % feature vector
            fv_train.eeg = squeeze(logvar_train.x);
            fv_test.eeg = squeeze(logvar_test.x);
            
            y_train  = group.ment(:,train);
            y_test   = vec2ind(group.ment(:,test));
            
            % train classifier
            C.eeg = train_RLDAshrink(fv_train.eeg,y_train); % shrinkage linear discriminant analysis
            
            % classification
            grouphat.eeg(foldIdx,:) = LDAmapping(C.eeg,fv_test.eeg); % use custom function
            cmat.eeg(:,:,foldIdx) = confusionmat(y_test, grouphat.eeg(foldIdx,:));
        end
        
        acc.eeg.ment(shiftIdx,stepIdx)   = trace((sum(cmat.eeg,3))) / sum(sum(sum(cmat.eeg,3),2),1);
    end
end

mean_acc.eeg.ment = mean(acc.eeg.ment,1)';

% display
time = (ival(:,2)/1000)';
figure('Name',subdir_list{vp}, 'Number','off')
plot(time, mean_acc.eeg.imag, 'b', time, mean_acc.eeg.ment,'r');
legend('motor imagery','mental arithmetic');
xlim([-5 17]); ylim([0.4 1]); grid on;
