%% Use unfold on the data
% used
% https://gin.g-node.org/agert/WildLab/src/master/code/unfold/unfold_freeviewing.m
% for input
clear variables
close all;
clc;
%% 
addpath('.../Matlab-resources/eeglab2020_0');
% load EEGlab
eeglab;

addpath('.../Matlab-resources/unfold-main') 
init_unfold
addpath('.../Matlab-resources/unfold-main/gramm')


basepath='.../processedData/village/preprocessed/';
cd(basepath);
cd('.../repos/wd-pilot-pipeline');
rec_vill = readtable('recordings_village.csv');

% save the data to use in Python
py_path = '.../data/village/processed/EEG';

%%
subs_to_include = [1,4,5,11,16,18,19,20,21,29,30,33,34,36,37,38,41,...
    42,43,44,45,46,47,48,49,50,53,54,56,57,58,59,60]; 

%%
cnt = 1;
bgrd = nan(length(subs_to_include),64,750);
body = nan(length(subs_to_include),64,750);
head = nan(length(subs_to_include),64,750);

for sub = 1:length(subjects)
    s = subjects(sub); 

    cd('.../recordedData/wd_village/');
    uidname = rec_vill{sub,1};
    uidname = uidname{1,1};

    savedata = [basepath, uidname, '/'];

    % add this new folder to the savedata path, so where the intermediate steps
    % will be saved
    savedata = [savedata, 'automated_preproc_new/'];
    %% Load the data
    if ismember(s,subs_to_include)
        cd(savedata);
        %% load & prepare the EEG file (do nto reject noisy segments)
        EEG = pop_loadset(sprintf('2a_cleanDataChannels_woRejection_%s.set',uidname),fullfile(savedata));

        % update sacc_amp/ sacc_dur
        % trg_path = '/net/store/nbp/projects/wd_ride_village/repos/wd-pilot-pipeline/data/village/processed/Trigger_MAD_newfile/';
        % %load(fullfile(trg_path,sprintf('sacc_amp_updated_%s.mat',uidname)));
        % load(fullfile(trg_path,sprintf('sacc_dur_%s.mat',uidname)));
        % 
        % for j = 1:length(EEG.event)
        %     EEG.event(j).saccade_amp = data(j);
        % end
        

        % redo the ica weights
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
        
        % redo interpolation 
        EEG_chan = pop_loadset(sprintf('1a_triggersFiltering_%s.set',uidname),fullfile(savedata));
        full_chanlocs = EEG_chan.chanlocs; % used for data cleaning and interpolation
        clear EEG_chan

        EEG = pop_interp(EEG,full_chanlocs,'spherical');
        EEG = eeg_checkset(EEG);
        
        %% prepare the event file
        
        % remove NaNs (for now)
        empty_events = false(size(EEG.event));
        for j = 1:length(EEG.event)
            event_fileds = fieldnames(EEG.event(1));
            for k = 1:numel(event_fileds)
                if isempty(EEG.event(j).(event_fileds{k}))
                    empty_events(j) = true;
                    break
                % for sacc_amp: detect empty events
                elseif isnan(EEG.event(j).(event_fileds{k}))
                    empty_events(j) = true;
                    break
                end
                
            end
        end
                
        EEG.event(empty_events) = [];
        
        % remove empty events == 3 (for now)        
        empty_events = false(size(EEG.event));
        for j = 1:length(EEG.event)
            if EEG.event(j).type == 3
                empty_events(j) = true;
            end
        end
        EEG.event(empty_events) = [];
        
        % add currentType and change type to fixation
        for j = length(EEG.event):-1:1
            EEG.event(j).currentType = EEG.event(j).type;
            EEG.event(j).type = 'fixation';
        end

        % change the numbers to meaningful strings
        for j = 1:length(EEG.event)
            switch EEG.event(j).currentType
                case 3 % atm not part of the events anymore
                    EEG.event(j).currentType = 'noName';
                case 2
                    EEG.event(j).currentType = 'background';
                case 1
                    EEG.event(j).currentType = 'body';
                case 0
                    EEG.event(j).currentType = 'head';
            end
        end
        
        % turn sacc_amp into float
        for j = 1:length(EEG.event)
            EEG.event(j).saccade_amp = str2double(EEG.event(j).saccade_amp);
        end
        
        EEG1 = EEG;
        %% Create design matrix
        cfgDesign               = [];
        % all events we have are fixations:
        cfgDesign.eventtypes    = {'fixation'}; 
        % since we do not have a clear group to compare the rest too: 
        cfgDesign.codingschema  = 'effects'; 
        % splines at quantiles of variable --> lot of data == lot of splines 
        %cfgDesign.splinespacing = 'quantile'; 
        % we have one categroical predictor and a non-linear predictor for
        % the saccade amplitude:
        %cfgDesign.formula       = 'y~1 + cat(currentType) + spl(saccade_amp,5)'; 
        cfgDesign.formula       = 'y~1 + cat(currentType)'; 
        % order them to have bgrd last (bgrd is reference)
        cfgDesign.categorical   = {'currentType',{'body','head','background'};};
        EEG                    = uf_designmat(EEG,cfgDesign);
        
        %% Timeshift
        cfgTimeshift            = [];
        cfgTimeshift.timelimits = [-0.5 1.0];
        EEG                    = uf_timeexpandDesignmat(EEG,cfgTimeshift);
        
        %% Reject noisy segments
        % we only do this, if we have noisy segments to reject
        if isfile(fullfile(savedata,sprintf('removed_intervals_%s.mat',uidname)))
            load(fullfile(savedata,sprintf('removed_intervals_%s.mat',uidname)));
            EEG                = uf_continuousArtifactExclude(EEG,'winrej',tmprej);
        end
        
        %% Fit the model
        % fit the model:
        cfgFit                  = [];
        cfgFit.precondition     = 1;
        cfgFit.lsmriterations   = 1500; % steps iterative solver should reach
        cfgFit.channel          = 1:length(EEG.chanlocs); % all channels
        EEG                     = uf_glmfit(EEG, cfgFit); 
        
        %% Make a massive uni-variate fit without de-convolution (Gert et al., 2022)
        if isfile(fullfile(savedata,sprintf('removed_intervals_%s.mat',uidname)))
            load(fullfile(savedata,sprintf('removed_intervals_%s.mat',uidname)));
            EEGepoch = uf_epoch(EEG,'winrej',tmprej,'timelimits',cfgTimeshift.timelimits);
        else
            EEGepoch = uf_epoch(EEG,'timelimits',cfgTimeshift.timelimits);
        end
        
        EEGepoch = uf_glmfit_nodc(EEGepoch);
        %% Get the betas
        % results condensed in new structure
        ufresultEp              = uf_condense(EEGepoch);
        ufresultEp.beta = ufresultEp.beta * 2; % * 2 due to effect coding with -1 / 1.
        
        ufresultEp = uf_predictContinuous(ufresultEp); % only overlap
        %ufresultEp = uf_predictContinuous(ufresultEp,'predictAt',{{'sac_amplitude',[0.5 1 2.5 5 7.5 10]},}); 
        ufresultEp = uf_addmarginal(ufresultEp);


        paramNames={ufresultEp.param.name};
        [~,paramPos_bg]=find(ismember(paramNames,'(Intercept)'));
        [~,paramPos_bo]=find(ismember(paramNames,'currentType_body'));
        [~,paramPos_he]=find(ismember(paramNames,'currentType_head'));
        

        bgrd(cnt,:,:) = ufresultEp.beta(:,:,paramPos_bg);
        body(cnt,:,:)  = ufresultEp.beta(:,:,paramPos_bo);
        head(cnt,:,:)  = ufresultEp.beta(:,:,paramPos_he);

        time = ufresultEp.times;
        
        %% Save the data
        
        save(fullfile(py_path,sprintf('betas_unfold_saccOnset_%s.mat',uidname)),'bgrd','body','head','time'); 
        cnt = cnt+1; % to avoid having nans in the data
    end
end
% save the data across all subjects
save(fullfile(py_path,'betasUnfold_saccOnset_allEl_justOverlap.mat'),...
     'bgrd','body','head','time'); 

chanLocs = EEG.chanlocs;
save(".../chanLocs.mat",'chanLocs');
