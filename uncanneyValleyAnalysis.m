%% Uncanny Valley analyses
%% Define paths
dataFolder = []; % Datapath to the csv file
%% Load the data
surveyResponses = readtable([dataFolder,filesep,'FacePerceptionResponses.csv']);
%% Extract demographics
age = surveyResponses.Age;
gender = surveyResponses.Gender;
familiarityVR = surveyResponses.AreYouFamiliarWithVREnvironmentsAndAvatars_;
%% Compute demographics
meanAge = mean(age);
stdAge = std(age);

females =  sum(contains(gender,'Female'));
males =  sum(contains(gender,'Male'));

veryFamiliar = sum(contains(familiarityVR,'1 - Very Familiar'));
familiar = sum(contains(familiarityVR,'2 - Familiar'));
neutral = sum(contains(familiarityVR,'3 - Neutral'));
unfamiliar = sum(contains(familiarityVR,'4 - Unfamiliar'));
veryUnfamiliar = sum(contains(familiarityVR,'5 - Very Unfamiliar'));
%% Determine order of questions
items = {'Artificial vs Natural','Human made vs Humanlike','No Definite Lifespam vs Mortal','Inanimate vs Livivng',...
    'Mechanical Movement vs Biological Movement','Synthetic vs Real','Reassuring vs Eerie','Numbring vs Freaky',...
    'Ordinary vs Supernatural','Bland vs Uncanny','Unemotional vs Hair-Raising','Uninspiring vs Spine-Tringling',...
    'Predictable vs Thrilling','Boring vs Shocking','Repulsive vs Agreeable','Ugly vs Beautiful','Messy vs Sleek',...
    'Crude vs Stylish','Very strange vs Very familiar'};
nItems = length(items);
%% Picture order
pictures = {'A10','A09','R10','R01','S02','S03','A08','S05','U07','A07','R06','U03','A04','R08','S10','R02',...
'R07','S09','U05','U01','R04','S04','U04','A01','A06','U10','U08','R03','U02','R09','S07','A05','U06','U09',...
'A02','A03','S08','S06','R05','S01'};
nPictures = length(pictures);
%%
clear surveyResponses
%% Load the data
surveyResponses = readmatrix([dataFolder,filesep,'FacePerceptionResponses.csv']);
%% Organize the data
ratingData = surveyResponses(:,5:end);
nSubjects = size(surveyResponses,1);
% Matrix Subjects x Questions x Picture
DATA = reshape(ratingData, nSubjects, nItems, nPictures);
% Remove participant 2 (was a test run)
DATA(2,:,:) = [];
%% Reorganize so all the groups of pictures are together
pictureType = cellfun(@(x) x(1), pictures, 'UniformOutput', false);  % Extracts the type (A, R, U, S)
pictureNumber = cellfun(@(x) str2double(x(2:end)), pictures); 
[~, sortOrder] = sortrows([string(pictureType)' pictureNumber']);

sortedPictureIDs = pictures(sortOrder);
sortedDATA = DATA(:, :, sortOrder);

% Compute the mean acrros subjects
meanOverSubs = squeeze(mean(sortedDATA,1));
%% Separate the different scales and do the plot of the mean per scale,
% in order to see better.
humannesIndex = 1:6;
eerinessIndex = 7:14;
attractivenessIndex = 15:18;
%% Prepare the data for analyses
humannesScale = sortedDATA(:,humannesIndex,:);
eerinessScale = sortedDATA(:,eerinessIndex,:);
attractivenesScale = sortedDATA(:,attractivenessIndex,:);
%% Statistical tests
% Humanness
meanHumanness = squeeze(mean(humannesScale(:,:,:),2));
reshaped_data = [squeeze(mean(meanHumanness(:,1:10),2)),squeeze(mean(meanHumanness(:,11:20),2)),...
    squeeze(mean(meanHumanness(:,21:30),2)),squeeze(mean(meanHumanness(:,31:end),2))]; 
% Perform the Friedman test
[p, tbl, stats] = friedman(reshaped_data, 1, 'off');
[c,m,h,gnames] = multcompare(stats);

% Plot the results
figure(1)
subplot(2,1,1)
boxplot(reshaped_data)
title('Humanness Across Participants')
set(gca,'XTickLabel',{'VR avatars','Realistics','Semi Realistics','Unrealistics'})
set(gca, 'FontSize', 14);
subplot(2,1,2)
multcompare(stats)
title('Humanness Differences')
yticklabels({'Unrealistics','Semi Realistics','Realistics','VR avatars'})
set(gca, 'FontSize', 14);

% Eeriness
meanEeriness = squeeze(mean(eerinessScale(:,:,:),2));
reshaped_data = [squeeze(mean(meanEeriness(:,1:10),2)),squeeze(mean(meanEeriness(:,11:20),2)),...
    squeeze(mean(meanEeriness(:,21:30),2)),squeeze(mean(meanEeriness(:,31:end),2))]; 
% Perform the Friedman test
[p, tbl, stats] = friedman(reshaped_data, 1, 'off');
[c,m,h,gnames] = multcompare(stats);

% Plot the results
figure(1)
subplot(2,1,1)
boxplot(reshaped_data)
title('Eeriness Across Participants')
set(gca,'XTickLabel',{'VR avatars','Realistics','Semi Realistics','Unrealistics'})
set(gca, 'FontSize', 14);
subplot(2,1,2)
multcompare(stats)
title('Eeriness Differences')
yticklabels({'Unrealistics','Semi Realistics','Realistics','VR avatars'})
set(gca, 'FontSize', 14);

% Attractiveness
meanAttractiveness = squeeze(mean(attractivenesScale(:,:,:),2));
reshaped_data = [squeeze(mean(meanAttractiveness(:,1:10),2)),squeeze(mean(meanAttractiveness(:,11:20),2)),...
    squeeze(mean(meanAttractiveness(:,21:30),2)),squeeze(mean(meanAttractiveness(:,31:end),2))]; 
% Perform the Friedman test
[p, tbl, stats] = friedman(reshaped_data, 1, 'off');
[c,m,h,gnames] = multcompare(stats);

% Plot the results
figure(1)
subplot(2,1,1)
boxplot(reshaped_data)
title('Attractiveness Across Participants')
set(gca,'XTickLabel',{'VR avatars','Realistics','Semi Realistics','Unrealistics'})
set(gca, 'FontSize', 14);
subplot(2,1,2)
multcompare(stats)
title('Attractiveness Differences')
yticklabels({'Unrealistics','Semi Realistics','Realistics','VR avatars'})
set(gca, 'FontSize', 14);
%% Plot figure for article
% First, average over the questions for each participant and picture:
avg_humanness = squeeze(mean(humannesScale, 2));
avg_eeriness  = squeeze(mean(eerinessScale, 2));
% Average over the 12 participants to obtain one value per picture:
pictureHumanness = mean(avg_humanness, 1); 
pictureEeriness  = mean(avg_eeriness, 1);
% Convert row vectors to column vectors:
pictureHumanness = pictureHumanness(:); 
pictureEeriness  = pictureEeriness(:);
% Prepare table for plotting
num_pictures = length(pictureHumanness);
group_vector = [ones(10,1); 2*ones(10,1); 3*ones(10,1); 4*ones(10,1)];
% Create a table with one row per picture
T_picture = table((1:num_pictures)', pictureHumanness, pictureEeriness, group_vector, ...
    'VariableNames', {'Picture','Humanness','Eeriness','Group'});
% Convert Group to categorical
T_picture.Group = categorical(T_picture.Group);
% Do the plot
figure;
hold on;
colors = [125, 23, 0;2, 49, 152;93, 192, 211;200, 180, 85]/255;
groupLabels = categories(T_picture.Group);

% Loop over each group to plot points with different colors
for i = 1:length(groupLabels)
    % Logical index for current group
    idx = T_picture.Group == groupLabels{i};
    
    % Scatter plot for current group with different colors
    scatter(T_picture.Humanness(idx), T_picture.Eeriness(idx), 100, colors(i,:), 'filled');

    scatter(mean(T_picture.Humanness(idx)), mean(T_picture.Eeriness(idx)), 400, colors(i,:),'*','LineWidth',3);
end
xlabel('Aggregated Humanness');
ylabel('Aggregated Eeriness');
legend('VR avatars','','Realistics','','Semi Realistic','','Unrealistics','', 'Location','Best');
grid on;
set(gca, 'FontSize', 16);

print(gcf, 'PaperFigureUncannyValley.png', '-dpng', '-r400')

