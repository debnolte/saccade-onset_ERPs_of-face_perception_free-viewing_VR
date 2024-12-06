%% TFCE ANOVA
%% 
clear variables
close all;
clc;

addpath('.../Matlab-resources/unfold-main') 
addpath('.../Matlab-resources/unfold-main/gramm')
addpath('.../Matlab-resources/unfold-main/lib/ept_TFCE/TFCE')
addpath('.../Matlab-resources/unfold-main/lib/ept_TFCE/TFCE/Dependencies')

%% get the data
cd ('.../TFCE_resources')
% where to save

py_path = '.../data/village/processed/EEG';

load(fullfile(py_path,sprintf('betasUnfold_saccOnset_allEl_justOverlap.mat')));
times = time;

thresh = 0.05;

Data{1,1} = bgrd(:,:,:);
Data{1,2} = body(:,:,:);  
Data{1,3} = head(:,:,:);
Data{2,1} = bgrd(:,:,:);
Data{2,2} = body(:,:,:); 
Data{2,3} = head(:,:,:);
%% 
% Set random stream depending on the clock value (unpredictable).
myRand = RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setGlobalStream(myRand);

nP = cell2mat(cellfun(@(x) size(x,1), Data, 'UniformOutput', false)); % Find the number of participants in each dataset

display ('Checking Data...')
if numel(unique(nP)) > 1
    error('All datsets must have the same number of participants (Balanced Design)');
end

nP   = unique(nP);
nCh  = size(Data{1},2);
nS   = size(Data{1},3);

% Load Channlocs
load(".../TFCE_resources/chanLocs.mat")
% Error Checking
if ~isequal(nCh, length(chanLocs))
    error ('Number of channels in data does not equal that of locations file')
end
display ('Done.')
%% Calculate the channels neighbours... using the modified version ChN2
display ('Calculating Channel Neighbours...')
ChN = ept_ChN2(chanLocs);
display ('Done.')
%% Create all variables in loop at their maximum size to increase performance
nPerm = 10000;
maxTFCE.A  = zeros(nPerm,1);
maxTFCE.B  = zeros(nPerm,1);
maxTFCE.AB = zeros(nPerm,1);
%% Calculate the actual T values of all data
% Calculate different T-values for mixed and repeated measures ANOVA
display ('Calculating Observed Statistics...')

F_Obs = ept_rmANOVA(Data);

E_H = [0.666, 1];

TFCE_Obs.A  = ept_mex_TFCE2D(double(F_Obs.A),  ChN, E_H);
TFCE_Obs.B  = ept_mex_TFCE2D(double(F_Obs.B),  ChN, E_H);
TFCE_Obs.AB = ept_mex_TFCE2D(double(F_Obs.AB), ChN, E_H);
display ('Done.')
%% Calculating the T value and TFCE enhancement of each different permutation
display ('Calculating Permutations...')
for i   = 1:nPerm

    F_Perm = ept_rmANOVA(Data,1);

    TFCE_Perm.A  = ept_mex_TFCE2D(F_Perm.A,  ChN, E_H);
    TFCE_Perm.B  = ept_mex_TFCE2D(F_Perm.B,  ChN, E_H);
    TFCE_Perm.AB = ept_mex_TFCE2D(F_Perm.AB, ChN, E_H);

    maxTFCE.A(i)  = max(max(abs(TFCE_Perm.A)));       % stores the maximum absolute value
    maxTFCE.B(i)  = max(max(abs(TFCE_Perm.B)));
    maxTFCE.AB(i) = max(max(abs(TFCE_Perm.AB)));
end
display ('Done.')
%% Calculating the p value from the permutation distribution
display ('Calculating Final Statistics...')

% add observed maximum
edges.A  = [maxTFCE.A;  max(max(abs(TFCE_Obs.A)))];
edges.B  = [maxTFCE.B;  max(max(abs(TFCE_Obs.B)))];
edges.AB = [maxTFCE.AB; max(max(abs(TFCE_Obs.AB)))];

[~,bin.A]      = histc(abs(TFCE_Obs.A),sort(edges.A));
P_Values.A     = 1-bin.A./(nPerm+2);
[~,bin.B]      = histc(abs(TFCE_Obs.B),sort(edges.B));
P_Values.B     = 1-bin.B./(nPerm+2);
[~,bin.AB]     = histc(abs(TFCE_Obs.AB),sort(edges.AB));
P_Values.AB    = 1-bin.AB./(nPerm+2);

display ('Done.')

%%
% find all indicies of significant clusters
idx = find(P_Values.B < thresh);
% convert linear indices to row and column subscripts
[Ch, S] = ind2sub(size(P_Values.B), idx);
% extract the corresponding p-values
p_values_below_threshold = P_Values.B(idx);


[min_P, ID] = min(P_Values.B(:));
[ChanMax, SMax]      = ind2sub(size(P_Values.B),ID);
max_Obs      = F_Obs.B(ID);

if min_P < thresh
    display(['Peak significance found at channel ', chanLocs(ChanMax).labels, ' at time point ', num2str(times(SMax)), ': T(', num2str(size(Data{1},1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);
else
    display('No siginificant clusters found!')
end

%% Find All Clusters
% Look for new clusters
x = TFCE_Obs.B; 
x(P_Values.B>thresh) = 0; % Threshold the data at alpha
Cp = ept_ClusRes(x, ChN, 0.01); % Calculate Positive Clusters
xn = x; % copy x
xn(x>0)=0; % Only show negative x
xn = abs(xn); % Make negative x's positive
Cn = ept_ClusRes(xn, ChN, 0.01); % Calculate Negative Clusters
C = Cn-Cp; % Combine the two
assignin('base', 'current_clusters', C);
b = unique(C); % How many different clusters are there?
b(b==0)=[]; % Eliminate the 0 from being a unique cluster

if numel(b)==0
    display (['There are no clusters of significant data at the p = ' num2str(thresh) ' threshold']);
    return
else
    ClusRes= cell(size(b,1),10);
    for i = 1:size(b,1)
        x = TFCE_Obs.B;
        x(C~=b(i))      = 0; %
        idPeak          = find(abs(x)==max(abs(x(:))));
        [PeakC, PeakS]  = ind2sub(size(x),idPeak);
        idSize          = find(C== b(i)); % find the rows and columns that are significant
        [SizeC, SizeS]  = ind2sub(size(C),idSize);
        
        ClusRes{i,1}    = PeakC(1); % peak channel (just the first of many possible peak channels (but averaging may result in a channel in between two that is not significant)!
        ClusRes{i,2}    = PeakS(1);
        ClusRes{i,3}    = F_Obs.B(PeakC(1),PeakS(1));
        ClusRes{i,4}    = P_Values.B(PeakC(1),PeakS(1));
        ClusRes{i,5}    = numel(idSize);
        ClusRes{i,6}    = numel(unique(SizeC));
        ClusRes{i,7}    = numel(unique(SizeS));
        ClusRes{i,8}   = [num2str(min(SizeS)), ' - ', num2str(max(SizeS))];
        % electrode:
        ClusRes{i,9}   = chanLocs(PeakC(1)).labels;
        % time interval
        ClusRes{i,10}   = [num2str(times(min(SizeS))*1000), ':', num2str(times(max(SizeS))*1000),'ms'];
    end
    [ClusRes{12,:}] = deal('idxChan','idxTime','Fval','pval','nrRowsColSig','nrChan','nrTime','timeIntlIdx','ChanName','timeInterval');
end


%% Save the Data
Info.Parameters.E_H         = E_H;
Info.Parameters.nPerm       = nPerm;
Info.Parameters.rSample     = 500;
Info.Parameters.type        = 'm';
Info.Parameters.nChannels   = nCh;
Info.Parameters.nSamples    = nS;
Info.Parameters.GroupSizes  = [size(Data{1,1}, 1); size(Data{2,1}, 1)];
Info.Electrodes.chanloc     = chanLocs;
Info.Electrodes.ChaNeigh    = ChN;

Results.Obs                 = F_Obs.B;
Results.TFCE_Obs            = TFCE_Obs;
Results.maxTFCE             = maxTFCE.B;
Results.P_Values            = P_Values.B;
Results.TFCE_Perm           = TFCE_Perm.B;
Results.idx                 = idx;
Results.Ch                  = Ch;
Results.S                   = S;
Results.p_vals_u_thresh     = p_values_below_threshold;

tfce_path = '.../TFCE_resources';
save(fullfile(tfce_path,'TFCEResults_Unfold_10000.mat'),'Info','Results','Data','ClusRes');

