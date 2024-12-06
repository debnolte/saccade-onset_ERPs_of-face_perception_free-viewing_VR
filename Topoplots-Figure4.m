%% Plot Topoplots without unfold
clear variables
close all;
clc;
%%
addpath('.../Matlab-resources/eeglab2020_0');
% load EEGlab
eeglab;
addpath('.../Matlab-resources/fieldtrip-20230215/');
ft_defaults

basepath='.../processedData/village/preprocessed/';
cd(basepath);
cd('.../repos/wd-pilot-pipeline');
rec_vill = readtable('recordings_village.csv');

face_path = '.../repos/wd-pilot-pipeline/FacePaper/images';

%%
subjects = [1,2,4,5,7,8,10,11,12,15,16,17,18,19,20,21,22,24,26,27,29,30,31,32,33,...
    34,36,37,38,41,42,43,44,45,46,47,48,49,50,51,53,54,55,56,57,58,59,60];

subs_to_include = [1,4,5,11,16,17,18,19,20,21,29,30,33,34,36,37,38,41,...
    42,43,44,45,46,47,48,49,50,51,53,54,56,57,58,59,60];

%% Load the data, saved in Topoplots_acrossSubj.m
% Load Channlocs
load(".../chanLocs.mat")

py_path = '.../data/village/processed/EEG';
load(fullfile(py_path,sprintf('FacePaper_ERPData_Full.mat')));
clear('gaze_bgrd','gaze_body','gaze_head','sacc_bgrd','sacc_body','sacc_head',...
    'diff_bgrd_body','diff_bgrd_head','diff_body_head');

avg_gaze = mean(avg_gaze,3);
avg_sacc = mean(avg_sacc,3);

% colormap:
cmi = [5, 48, 97; 30, 97, 164 ;60, 138, 190; 123, 182, 214; 187, 218, 233;...
    230, 239, 244; 250, 234, 225; 250, 200, 174; 235, 145, 114; 207, 82, 70;...
    171, 22, 42; 103, 0, 31];

%% Plot Gaze onset topoplots
lg = [[-60 -40];[90 110];[170 190]];

for i = 1:3
    gaze = mean(avg_gaze(:,find(times>=lg(i,1),1):find(times<=lg(i,2),1,'last')),2);

    figure;
    topoplot(gaze,chanLocs,                       ...% plot TFCE_Obs or pPlot
        'electrodes','on',                              ...% display markers ("labels" shows electrode names
        'whitebk','on',                                     ...
        'colormap',cmi,                                     ...
        'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
        'style','map',                                      ...
        'numcontour',6,                                     ...% Increase for more contour details
        'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
        'maplimits',[-2,2]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
    set(gcf, 'PaperOrientation', 'landscape')
    titleStr = strcat('Topoplot_gaze_',string(lg(i,1)),'_',string(lg(i,2)),'.png');
    saveas(gcf,fullfile(face_path,titleStr))
end

%% Plot Saccade onset topoplots
lg = [[0 20];[150 170];[230 250]];

for i = 1:3
    sacc = mean(avg_sacc(:,find(times>=lg(i,1),1):find(times<=lg(i,2),1,'last')),2);

    figure;
    topoplot(sacc,chanLocs,                       ...% plot TFCE_Obs or pPlot
        'electrodes','on',                              ...% display markers ("labels" shows electrode names
        'whitebk','on',                                     ...
        'colormap',cmi,                                     ...
        'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
        'style','map',                                      ...
        'numcontour',6,                                     ...% Increase for more contour details
        'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
        'maplimits',[-2,2]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
    set(gcf, 'PaperOrientation', 'landscape')
    titleStr = strcat('Topoplot_sacc_',string(lg(i,1)),'_',string(lg(i,2)),'.png');
    saveas(gcf,fullfile(face_path,titleStr))
end

%% Plot the Colorbar
figure;
topoplot(gaze,chanLocs,                       ...% plot TFCE_Obs or pPlot
    'electrodes','on',                              ...% display markers ("labels" shows electrode names
    'whitebk','on',                                     ...
    'colormap',cmi,                                     ...
    'shading','interp',                                 ...% (flat or interp) useless with style is "fill"
    'style','map',                                      ...
    'numcontour',6,                                     ...% Increase for more contour details
    'plotrad',max(abs(cell2mat({chanLocs.radius}))),    ...
    'maplimits',[-2,2]);                                   % 'map' 'contour' 'both' 'fill' 'blank'
c = colorbar;
c.Limits = ([-2,2]);
c.Ticks = [-2 -1 0 1 2];
c.FontSize = 14;
ylabel(c,'\muV','FontSize',18,'Rotation',90)
set(gca,'fontname','arial')
set(gcf, 'PaperOrientation', 'landscape')
titleStr = strcat('Colorbar_AvgPlots.png');
saveas(gcf,fullfile(face_path,titleStr))

