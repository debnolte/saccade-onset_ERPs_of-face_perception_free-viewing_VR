%% Collect the data of all participants to be unsed for the final plots
clear variables
close all;
clc;
%%
addpath('.../Matlab-resources/eeglab2020_0');
% load EEGlab
eeglab;

basepath='.../processedData/village/preprocessed/';
cd(basepath);
cd('.../repos/wd-pilot-pipeline');
rec_vill = readtable('recordings_village.csv');
%%
subs_to_include = [1,4,5,11,16,18,19,20,21,29,30,33,34,36,37,38,41,...
    42,43,44,45,46,47,48,49,50,53,54,56,57,58,59,60]; 

%% Create data to save
% gaze onset:
avg_gaze  = nan(64,400,length(subs_to_include)); % average across all gaze trials
gaze_bgrd = nan(64,400,length(subs_to_include)); % bgrd
gaze_body = nan(64,400,length(subs_to_include)); % body
gaze_head = nan(64,400,length(subs_to_include)); % head
% saccade onset:
avg_sacc  = nan(64,400,length(subs_to_include)); % average across all saccade trials
sacc_bgrd = nan(64,400,length(subs_to_include)); % bgrd
sacc_body = nan(64,400,length(subs_to_include)); % body
sacc_head = nan(64,400,length(subs_to_include)); % head
% differences for saccade consitions:
diff_bgrd_body = nan(64,400,length(subs_to_include)); % bgrd - body
diff_bgrd_head = nan(64,400,length(subs_to_include)); % bgrd - head
diff_body_head = nan(64,400,length(subs_to_include)); % body - head

%% loop through all subjects & collect the data
cnt = 1; % counted up to save without nans in between
for sub = 1:length(subjects)
    s = subjects(sub); 
    cd('.../recordedData/wd_village/');
    uidname = rec_vill{sub,1};
    uidname = uidname{1,1};
    savedata = [basepath, uidname, '/'];
    savedata = [savedata, 'automated_preproc_new/'];
    %% Load the data
    if ismember(s,subs_to_include)
        cd(savedata);

        %% gaze
        EEG = pop_loadset(sprintf('4acorr_gaze_interpolation_%s.set',uidname),fullfile(savedata));
        % epoch for 3 conditions: 2 = 'bgrd'; 1 = 'body'; 0 = 'face'
        EEG_bgrd = pop_epoch(EEG, {'2'}, [-0.3 0.5]);
        EEG_body = pop_epoch(EEG, {'1'}, [-0.3 0.5]);
        EEG_head = pop_epoch(EEG, {'0'}, [-0.3 0.5]);
        % epoch across all trials
        EEG = pop_epoch(EEG, {}, [-0.3 0.5]);
        
        % add it to the matrices
        avg_gaze(:,:,cnt)  = mean(EEG.data,3); % all
        gaze_bgrd(:,:,cnt) = mean(EEG_bgrd.data,3); % bgrd
        gaze_body(:,:,cnt) = mean(EEG_body.data,3); % body
        gaze_head(:,:,cnt) = mean(EEG_head.data,3); % head

        %% saccade
        EEG = pop_loadset(sprintf('4a_interpolation_%s.set',uidname),fullfile(savedata));
        % epoch for 3 conditions: 2 = 'bgrd'; 1 = 'body'; 0 = 'face'
        EEG_bgrd = pop_epoch(EEG, {'2'}, [-0.3 0.5]);
        EEG_body = pop_epoch(EEG, {'1'}, [-0.3 0.5]);
        EEG_head = pop_epoch(EEG, {'0'}, [-0.3 0.5]);
        % epoch across all trials
        EEG = pop_epoch(EEG, {}, [-0.3 0.5]);
        
        % add it to the matrices
        avg_sacc(:,:,cnt)  = mean(EEG.data,3); % all
        sacc_bgrd(:,:,cnt) = mean(EEG_bgrd.data,3); % bgrd
        sacc_body(:,:,cnt) = mean(EEG_body.data,3); % body
        sacc_head(:,:,cnt) = mean(EEG_head.data,3); % head

        %% differences for saccade conditions 
        diff_bgrd_body(:,:,cnt) = mean(EEG_bgrd.data,3) - mean(EEG_body.data,3); % bgrd-body
        diff_bgrd_head(:,:,cnt) = mean(EEG_bgrd.data,3) - mean(EEG_head.data,3); % bgrd-head
        diff_body_head(:,:,cnt) = mean(EEG_body.data,3) - mean(EEG_head.data,3); % body-head
        %% count cnt up
        cnt = cnt + 1;
    end 
end
%% save the data
py_path = '...data/village/processed/EEG';
times = EEG_bgrd.times; % is the same for all subjctes so we only need it once
save(fullfile(py_path,sprintf('FacePaper_ERPData_Full.mat')),'avg_gaze','gaze_bgrd',...
    'gaze_body','gaze_head','avg_sacc','sacc_bgrd','sacc_body','sacc_head',...
    'diff_bgrd_body','diff_bgrd_head','diff_body_head','times');



%% Save the individual subject data
currElec = 'PO7';
el_idx = find(strcmp({EEG.chanlocs.labels}, currElec) == 1); %find the index of electrode

sub = 42; % subject to run (can be replaced by for loop for all subjects)
s = subjects(sub); 
cd('.../recordedData/wd_village/');
uidname = rec_vill{sub,1};
uidname = uidname{1,1};
savedata = [basepath, uidname, '/'];
savedata = [savedata, 'automated_preproc_new/'];
cd(savedata);

EEG = pop_loadset(sprintf('4acorr_gaze_interpolation_%s.set',uidname),fullfile(savedata));
% get the gaze avg for each condition
EEG_bgrd = pop_epoch(EEG, {'2'}, [-0.3 0.5]);
EEG_body = pop_epoch(EEG, {'1'}, [-0.3 0.5]);
EEG_head = pop_epoch(EEG, {'0'}, [-0.3 0.5]);
gaze_bgrd = mean(EEG_bgrd.data(el_idx,:,:), 3);
gaze_body = mean(EEG_body.data(el_idx,:,:), 3);
gaze_head = mean(EEG_head.data(el_idx,:,:), 3);

% get data for trial plot
EEG = pop_epoch(EEG, {}, [-0.3 0.5]);
data = squeeze(EEG.data(el_idx,:,:));
times = EEG.times;
event = EEG.event;

% get the sacc avg for each condition
EEG = pop_loadset(sprintf('4a_interpolation_%s.set',uidname),fullfile(savedata));
EEG_bgrd = pop_epoch(EEG, {'2'}, [-0.3 0.5]);
EEG_body = pop_epoch(EEG, {'1'}, [-0.3 0.5]);
EEG_head = pop_epoch(EEG, {'0'}, [-0.3 0.5]);
sacc_bgrd = mean(EEG_bgrd.data(el_idx,:,:), 3);
sacc_body = mean(EEG_body.data(el_idx,:,:), 3);
sacc_head = mean(EEG_head.data(el_idx,:,:), 3);

% save everything
py_path = '.../data/village/processed/EEG';
save(fullfile(py_path,sprintf('Sub%s_trialOrderSaccAmp.mat',string(sub))),...
    'data','times','event','gaze_bgrd','gaze_body','gaze_head',"sacc_bgrd","sacc_body","sacc_head");

