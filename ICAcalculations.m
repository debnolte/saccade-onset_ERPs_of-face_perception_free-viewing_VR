%% Set paths
dataFolder = []; % Datapath to folder with ICA results.
%% Get participants
tmp = dir(fullfile(dataFolder));
participants = [];
inx = 1;

for pId = 1:size(tmp,1)
    if tmp(pId).name(1) == '.' || contains(tmp(pId).name,'.set') || contains(tmp(pId).name,'.fdt')
        continue
    else
        participants(inx).name = tmp(pId).name;
        participants(inx).folder = tmp(pId).folder;
        participants(inx).date = tmp(pId).date;
        inx = inx + 1;
    end
end
clear tmp
%% 
for sub = 1:length(participants)
    EEGOUT = pop_loadset(sprintf('2a_cleanDataChannels_%s.set',participants(sub).name),fullfile(dataFolder));
    % Load the ICA outputs
    ICAfolder = [dataFolder,filesep,participants(sub).name];
    mod = loadmodout15(ICAfolder);
    % Apply ICA weights to data
    EEGOUT.icasphere = mod.S;
    EEGOUT.icaweights = mod.W;
    EEGOUT.icawinv = [];
    EEGOUT.icaact = [];
    EEGOUT.icachansind = [];
    EEGOUT = eeg_checkset(EEGOUT);
    % Use iclabel to determine which ICs to reject
    EEGOUT = iclabel(EEGOUT);
    
    % List components that should be rejected
    components2remove = [];
    muscle = [];
    eye = [];
    heart = [];
    line = [];
    channel = [];
    for component = 1:length(EEGOUT.chanlocs)-1
         % Muscle
         if EEGOUT.etc.ic_classification.ICLabel.classifications(component,2) > .80
             components2remove = [components2remove component];
             muscle = [muscle component];
         end
         % Eye
         if EEGOUT.etc.ic_classification.ICLabel.classifications(component,3) > .9
             components2remove = [components2remove component];
             eye = [eye component];
         end
         % Heart
         if EEGOUT.etc.ic_classification.ICLabel.classifications(component,4) > .9
             components2remove = [components2remove component];
             heart = [heart component];
         end
         % Line noise
         if EEGOUT.etc.ic_classification.ICLabel.classifications(component,5) > .9
             components2remove = [components2remove component];
             line = [line component];
         end
         % Channel noise
         if EEGOUT.etc.ic_classification.ICLabel.classifications(component,6) > .9
             components2remove = [components2remove component];
             channel = [channel component];
         end
    end
    allSubsMuscle(sub) = {muscle};
    allSubsEye(sub) = {eye};
    allSubsHeart(sub) = {heart};
    allSubsLine(sub) = {line};
    allSubsChannel(sub) = {channel};
end

for i = 1:length(allSubsChannel)
    chanComponents(i) = length(allSubsChannel{i});
    lineComponents(i) = length(allSubsLine{i});
    eyeComponents(i) = length(allSubsEye{i});
    heartComponents(i) = length(allSubsHeart{i});
    muscleComponents(i) = length(allSubsMuscle{i});
end

meanChan = mean(chanComponents);
stdChan = std(chanComponents);

meanLine = mean(lineComponents);
stdLine = std(lineComponents);

meanEye = mean(eyeComponents);
stdEye = std(eyeComponents);

meanHeart = mean(heartComponents);
stdHeart = std(heartComponents);

meanMuscle = mean(muscleComponents);
stdMuscle = std(muscleComponents);

% Define data
Names = {'Channel Noise'; 'Line Noise'; 'Eye'; 'Heart';'Muscle'};
mean = [meanChan;meanLine;meanEye;meanHeart;meanMuscle];
std = [stdChan;stdLine;stdEye;stdHeart;stdMuscle];
% Create a table
T = table(Names, mean, std);
% Display the table
disp(T);
