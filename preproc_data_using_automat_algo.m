function [EEG] = preproc_data_using_automat_algo(uidname,trgname, filt)
%% This function loads the EEG data and applies all preprocessing steps

% Inputs
% uidname   - subject uid, string
% trgname   - triggerfile name
% filt      - bool, True if you want to filter again, False otherwise
% filt_val  - value you want to filter with (for high-pass)
%
% Outputs
% EEG       - preprocessed EEG data, in eeglab form
%% Load necessary code
basepath='/net/store/nbp/projects/wd_ride_village/processedData/village/preprocessed/';
cd(basepath);

savedata = [basepath, uidname, '/'];

% add this new folder to the savedata path, so where the intermediate steps
% will be saved
savedata = [savedata, 'automated_preproc_new/'];
cd(savedata);

%trgname = 'TRF_gazeOnset_pos_shift_';
%% Load data
EEG     = pop_loadset(sprintf('2a_cleanDataChannels_woRejection_%s.set',uidname),fullfile(savedata));
EEG     = eeg_checkset(EEG);

%% Reload trigger file
trgpath='/net/store/nbp/projects/wd_ride_village/repos/wd-pilot-pipeline/data/village/processed/Trigger_MAD_newfile';
%trgpath='/net/store/nbp/projects/wd_ride_village/repos/wd-pilot-pipeline/data/village/processed/TriggerInfo_MAD';
EEG = pop_importevent(EEG,'event',fullfile(trgpath,strcat(trgname,uidname,'.csv')),'fields', {'latency', 'type'}, 'skipline', 1,'append','no');

%EEG = pop_importevent(EEG,'event',fullfile(trgpath,strcat(trgname,uidname,'.csv')),'fields', {'latency', 'type'}, 'skipline', 1,'append','no');
EEG=eeg_checkset(EEG);

%% Exclude data segements
if isfile(fullfile(savedata,sprintf('removed_intervals_%s.mat',uidname)))
    load(fullfile(savedata,sprintf('removed_intervals_%s.mat',uidname)));
    EEG = eeg_eegrej(EEG,tmprej);
end
EEG = eeg_checkset(EEG);

%% Exclude ICs
load(fullfile(savedata,sprintf('removed_components_%s.mat',uidname)));

outDir = fullfile(savedata, 'amica');
mod = loadmodout15(outDir);

EEG.icaweights  = mod.W;
EEG.icasphere   = mod.S;
EEG.icawinv     = [];
EEG.icaact      = [];
EEG.icachansind = [];
EEG             = eeg_checkset(EEG);
% get the bad components out of the saved ICA file and re-reject them 
EEG             = pop_subcomp(EEG,components_to_remove);

% save the results
EEG = eeg_checkset(EEG);

clear ic mod;

%% Interpolating missing channels
% all channels that are not needed
EEG_chan = pop_loadset(sprintf('1a_triggersFiltering_%s.set',uidname),fullfile(savedata));
full_chanlocs = EEG_chan.chanlocs; % used for data cleaning and interpolation
clear EEG_chan

EEG = pop_interp(EEG,full_chanlocs,'spherical');

EEG = eeg_checkset(EEG);
 % check if duplicate channel label
if isfield(EEG.chanlocs, 'labels')
    if length( { EEG.chanlocs.labels } ) > length( unique({ EEG.chanlocs.labels } ) )
        disp('Warning: some channels have the same label');
    end
end

% save the results
EEG = eeg_checkset(EEG);

clear tmp alldel idxs;
%% High Pass filter
EEG_filt = pop_eegfiltnew(EEG, filt, []);
%% Epoch the data and add filt to EEG
EEG = pop_epoch(EEG, {}, [-0.2 0.5]);
EEG_filt = pop_epoch(EEG_filt, {}, [-0.2 0.5]);
EEG.filt = EEG_filt.data;
%%
basepath='/net/store/nbp/projects/wd_ride_village/Analysis/FacePaper';
cd(basepath);

end