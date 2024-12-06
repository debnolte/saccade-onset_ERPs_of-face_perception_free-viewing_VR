%% 
clear variables
close all;
clc;

addpath('.../Matlab-resources/unfold-main') 
addpath('.../Matlab-resources/unfold-main/gramm')
addpath('.../Matlab-resources/unfold-main/lib/ept_TFCE/TFCE')
addpath('.../Matlab-resources/unfold-main/lib/ept_TFCE/TFCE/Dependencies')
addpath('.../Matlab-resources/eeglab2020_0');

%% Load the Data
tfce_path = '.../TFCE_resources';
load(fullfile(tfce_path,'TFCEResults_Unfold_10000.mat'));
% update tfce_path to the plotting folder
tfce_path = '.../TFCE_resources/TFCE_Plots';

% load the time vector
py_path = '.../data/village/processed/EEG';
load(fullfile(py_path,sprintf('betasUnfold_saccOnset_allEl_justOverlap.mat')));
times = time;
clear('py_path','time','bgrd','body','head')
times = times*1000; % update times for plotting

% Load Channlocs
load(".../TFCE_resources/chanLocs.mat")

%% Figure 5: Activation Maps
% Time points to plot:
timepoints = [20,106,150,200,240];

% Prepare plotting
bgrd_body = squeeze(mean(Data{1,1} - Data{1,2}));
bgrd_head = squeeze(mean(Data{1,1} - Data{1,3}));
body_head = squeeze(mean(Data{1,2} - Data{1,3})); 
diff      = {bgrd_body,bgrd_head,body_head};
titles    = {'bgrd_body','bgrd_head','body_head'};


cm = [0 0 1; 1 1 1; 1 0 0];                   % Basic Colormap
cmi       = interp1([-6; 0; 6], cm, (-6:6));  % interpolated Color

% loop through the timepoints
for j = 1:5
    idx = find(times == timepoints(j)); 

    % loop through the three conditions
    for i = 1:3
        diff_i = diff{i};
        figure;
        topoplot(diff_i(:,idx)',chanLocs,                       ...% plot TFCE_Obs or pPlot
            'electrodes','on',                              ...% display markers ("labels" shows electrode names
            'whitebk','on',                                     ...
            'colormap',cmi,                                     ...
            'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
            'style','map',                                      ...
            'numcontour',6,                                     ...% Increase for more contour details
            'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
            'maplimits',[-1,1]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
        set(gcf, 'PaperOrientation', 'landscape')
        saveas(gcf,fullfile(tfce_path,sprintf('Fig5_%s_%dms_topoDiff.png',titles{i},timepoints(j))))
        close(gcf);
    end
end

% Plot the topoplot for the colorbar
figure;
topoplot(diff_i(:,idx)',chanLocs,                       ...% plot TFCE_Obs or pPlot
    'electrodes','on',                              ...% display markers ("labels" shows electrode names
    'whitebk','on',                                     ...
    'colormap',cmi,                                     ...
    'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
    'style','map',                                      ...
    'numcontour',6,                                     ...% Increase for more contour details
    'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
    'maplimits',[-1,1]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
c = colorbar;
c.Limits = ([-1,1]);
c.Ticks = [-1 -0.5 0 0.5 1];
c.FontSize = 14;
set(gca,'fontname','arial')
set(gcf, 'PaperOrientation', 'landscape')
saveas(gcf,fullfile(tfce_path,'Fig5_Colorbar_Diffplots_TFCE.png'))

%% Figure 6: Plot the effect size of the cluster
bgrd_body = Data{1,1} - Data{1,2};
bgrd_head = Data{1,1} - Data{1,3};
body_head = Data{1,2} - Data{1,3}; 

% Shorten the arrays to contain the relevant time intervals only:
cluterStart = 260; % hardcoded, taken from ClusRes
clusterEnd  = 374; % hardcoded, taken from ClusRes
bgrd_body = bgrd_body(:,:,cluterStart:clusterEnd);
bgrd_head = bgrd_head(:,:,cluterStart:clusterEnd);
body_head = body_head(:,:,cluterStart:clusterEnd);
time = times(cluterStart:clusterEnd);


% Get the index of the absolute maximum value for each electrode and each
% participant 
[~, idx_bgrd_body] = max(abs(bgrd_body), [], 3);
[~, idx_bgrd_head] = max(abs(bgrd_head), [], 3);
[~, idx_body_head] = max(abs(body_head), [], 3);

idx_bgrd_body = sub2ind(size(bgrd_body), ...
    repmat((1:size(bgrd_body,1))',[1,size(bgrd_body,2)]),...
    repmat(1:size(bgrd_body,2),[size(bgrd_body,1),1]),...
    idx_bgrd_body);
idx_bgrd_head = sub2ind(size(bgrd_head), ...
    repmat((1:size(bgrd_head,1))',[1,size(bgrd_head,2)]),...
    repmat(1:size(bgrd_head,2),[size(bgrd_head,1),1]),...
    idx_bgrd_head);
idx_body_head = sub2ind(size(body_head), ...
    repmat((1:size(body_head,1))',[1,size(body_head,2)]),...
    repmat(1:size(body_head,2),[size(body_head,1),1]),...
    idx_body_head);


% Get the actual values associated with the generated indices
max_bgrd_body = bgrd_body(idx_bgrd_body);
max_bgrd_head = bgrd_head(idx_bgrd_head);
max_body_head = body_head(idx_body_head);

% Get the avgerage across participants for each electrode 
avg_bgrd_body = mean(max_bgrd_body,1);
avg_bgrd_head = mean(max_bgrd_head,1);
avg_body_head = mean(max_body_head,1);
avg_cond = {avg_bgrd_body,avg_bgrd_head,avg_body_head};

% Get the standard error of mean across participants for each electrode 
stderror_bgrd_body = std(max_bgrd_body,1)/sqrt(size(max_bgrd_body,1));
stderror_bgrd_head = std(max_bgrd_head,1)/sqrt(size(max_bgrd_head,1));
stderror_body_head = std(max_body_head,1)/sqrt(size(max_body_head,1));
stderror_cond = {stderror_bgrd_body,stderror_bgrd_head,stderror_body_head};

% Get the corresponding time points
[~, idx_bgrd_body] = max(abs(bgrd_body), [], 3);
[~, idx_bgrd_head] = max(abs(bgrd_head), [], 3);
[~, idx_body_head] = max(abs(body_head), [], 3);

time_bgrd_body = time(idx_bgrd_body);
time_bgrd_head = time(idx_bgrd_head);
time_body_head = time(idx_body_head);

% Average time over participants
avg_time_bgrd_body = mean(time_bgrd_body,1);
avg_time_bgrd_head = mean(time_bgrd_head,1);
avg_time_body_head = mean(time_body_head,1);
avg_time_cond = {avg_time_bgrd_body,avg_time_bgrd_head,avg_time_body_head};

% Standard error of mean of the time over participants
stderror_time_bgrd_body = std(time_bgrd_body,1)/sqrt(size(time_bgrd_body,1));
stderror_time_bgrd_head = std(time_bgrd_head,1)/sqrt(size(time_bgrd_head,1));
stderror_time_body_head = std(time_body_head,1)/sqrt(size(time_body_head,1));
stderror_time_cond = {stderror_time_bgrd_body,stderror_time_bgrd_head,stderror_time_body_head};

% Plot the results
titles    = {'bgrd_body','bgrd_head','body_head'};
cm = [0 0 1; 1 1 1; 1 0 0];                   % Basic Colormap
cmi       = interp1([-6; 0; 6], cm, (-6:6));  % interpolated Color
cm = [1 1 1; 0 0 1];                   % Basic Colormap
cmi1 = interp1([0; 6], cm, (0:6));  % interpolated Color

% loop through the three conditions
for i = 1:3
    % avg topoplot
    avg_i = avg_cond{i};
    figure;
    topoplot(avg_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
        'electrodes','on',                              ...% display markers ("labels" shows electrode names
        'whitebk','off',                                     ...
        'colormap',cmi,                                     ...
        'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
        'style','map',                                      ...
        'numcontour',6,                                     ...% Increase for more contour details
        'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
        'maplimits',[-1.8,1.8]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
    set(gcf, 'PaperOrientation', 'landscape')
    saveas(gcf,fullfile(tfce_path,sprintf('Fig6_%s_avgVals.png',titles{i})))
    close(gcf);
    
    % stderror topoplot
    std_i = stderror_cond{i};
    figure;
    topoplot(std_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
        'electrodes','on',                              ...% display markers ("labels" shows electrode names
        'whitebk','off',                                     ...
        'colormap',cmi1,                                     ...
        'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
        'style','map',                                      ...
        'numcontour',6,                                     ...% Increase for more contour details
        'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
        'maplimits',[0,0.6]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
    set(gcf, 'PaperOrientation', 'landscape')
    saveas(gcf,fullfile(tfce_path,sprintf('Fig6_%s_stderrorVals.png',titles{i})))
    close(gcf);

    % average time points
    time_i = avg_time_cond{i};
    figure;
    topoplot(time_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
        'electrodes','on',                              ...% display markers ("labels" shows electrode names
        'whitebk','off',                                     ...
        'colormap',cmi,                                     ...
        'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
        'style','map',                                      ...
        'numcontour',6,                                     ...% Increase for more contour details
        'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
        'maplimits',[18,246]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
    set(gcf, 'PaperOrientation', 'landscape')
    saveas(gcf,fullfile(tfce_path,sprintf('Fig6_%s_avgTime.png',titles{i})))
    close(gcf);

    % stderror time points
    s_time_i = stderror_time_cond{i};
    figure;
    topoplot(s_time_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
        'electrodes','on',                              ...% display markers ("labels" shows electrode names
        'whitebk','off',                                     ...
        'colormap',cmi1,                                     ...
        'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
        'style','map',                                      ...
        'numcontour',6,                                     ...% Increase for more contour details
        'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
        'maplimits',[7,15]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
    set(gcf, 'PaperOrientation', 'landscape')
    saveas(gcf,fullfile(tfce_path,sprintf('Fig6_%s_stderrorTime.png',titles{i})))
    close(gcf);
end

% Plot the topoplot for the colorbars
% avg topoplot
figure;
topoplot(avg_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
    'electrodes','on',                              ...% display markers ("labels" shows electrode names
    'whitebk','off',                                     ...
    'colormap',cmi,                                     ...
    'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
    'style','map',                                      ...
    'numcontour',6,                                     ...% Increase for more contour details
    'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
    'maplimits',[-1.8,1.8]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
c = colorbar;
c.Limits = ([-1.8,1.8]);
c.Ticks = [-1.8 -0.9 0 0.9 1.8];
c.FontSize = 14;
ylabel(c,'\muV','FontSize',18,'Rotation',90)
set(gca,'fontname','arial')
set(gcf, 'PaperOrientation', 'landscape')
saveas(gcf,fullfile(tfce_path,'Fig6_Colorbar_avgVals.png'))

% stderror topoplot
figure;
topoplot(std_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
    'electrodes','on',                              ...% display markers ("labels" shows electrode names
    'whitebk','off',                                     ...
    'colormap',cmi1,                                     ...
    'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
    'style','map',                                      ...
    'numcontour',6,                                     ...% Increase for more contour details
    'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
    'maplimits',[0,0.6]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
c = colorbar;
c.Limits = ([0,0.6]);
c.Ticks = [0 0.15 0.3 0.45 0.6];
c.FontSize = 14;
ylabel(c,'\muV','FontSize',18,'Rotation',90)
set(gca,'fontname','arial')
set(gcf, 'PaperOrientation', 'landscape')
saveas(gcf,fullfile(tfce_path,'Fig6_Colorbar_stderrorVals.png'))


% average time points
figure;
topoplot(time_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
    'electrodes','on',                              ...% display markers ("labels" shows electrode names
    'whitebk','off',                                     ...
    'colormap',cmi,                                     ...
    'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
    'style','map',                                      ...
    'numcontour',6,                                     ...% Increase for more contour details
    'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
    'maplimits',[18,246]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
c = colorbar;
c.Limits = ([18,246]);
c.Ticks = [18 75 132 189 246];
c.FontSize = 14;
ylabel(c,'time','FontSize',18,'Rotation',90)
set(gca,'fontname','arial')
set(gcf, 'PaperOrientation', 'landscape')
saveas(gcf,fullfile(tfce_path,'Fig6_Colorbar_avgTime.png'))

% stderror time points
figure;
topoplot(s_time_i(:)',chanLocs,                       ...% plot TFCE_Obs or pPlot
    'electrodes','on',                              ...% display markers ("labels" shows electrode names
    'whitebk','off',                                     ...
    'colormap',cmi1,                                     ...
    'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
    'style','map',                                      ...
    'numcontour',6,                                     ...% Increase for more contour details
    'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
    'maplimits',[7,15]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
c = colorbar;
c.Limits = ([7,15]);
c.Ticks = [7 9 11 13 15];
c.FontSize = 14;
ylabel(c,'time','FontSize',18,'Rotation',90)
set(gca,'fontname','arial')
set(gcf, 'PaperOrientation', 'landscape')
saveas(gcf,fullfile(tfce_path,'Fig6_Colorbar_stderrorTime.png'))

