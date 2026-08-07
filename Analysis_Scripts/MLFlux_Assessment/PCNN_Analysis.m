%% Load Data
model = 'Simple'; % Change this to Simple, UpperGly, Glycolysis, GlyPPP, or CCM

% Files can be reproduced via Data_Generation_and_Simulation
predictedLabels = load(strcat(model,'_FullPredictedLabels.txt'));
trueLabels = load(strcat('../../Test_Full_MLFlux/TestingData/',model,'_Labels_100_Seed0.txt'));
labelsToKeep = readmatrix(strcat('../../../../InformationSheets/',model,'_Labels_to_Keep.xlsx'));
labelsToKeep = find(labelsToKeep == 1);

%% Process labels to keep (an extra step if doing CCM 3Tracer)
indices = reshape(1:4800, 160, []);
indices = indices(1:24,:);
indices = indices(:);
labelsToKeep = labelsToKeep(indices);
labelsToKeep = find(labelsToKeep == 1);

%% Get sum of labels in each row
sumPredictedLabels = sum(predictedLabels,2);
trueSum = mean(sum(trueLabels,2));

figure(1)
clf
histogram(sumPredictedLabels)

xline(trueSum,'-',{'Expected Sum',trueSum})

%% Get the MAE of the labels
trueFullLabels = trueLabels(repmat(1:end,100,1),:);
mae_IndividaulLabels = mean(abs(trueFullLabels(:,labelsToKeep) - predictedLabels(:,labelsToKeep)));
mae_FullLabels = mean(abs(trueFullLabels(:,labelsToKeep) - predictedLabels(:,labelsToKeep)),2);

%% Plot the sum versus MAE
figure(3)
clf
scatterhist(sumPredictedLabels-trueSum,mae_FullLabels,'NBins',[20,20]);
set(gcf,'Position',[100,100,800,800]);
xlabel('Sum of predicted labels')
ylabel('Mean absolute error of labeling')
xlim([-4,6.5])
ylim([0,0.04])
set(gca,'TickDir','out')

%% Plot color density chart
figure(4)
clf
f=densityScatterChart(sumPredictedLabels-trueSum,mae_FullLabels);
set(gcf,'Position',[100,100,800,800]);
xlabel('Sum of predicted labels')
ylabel('Mean absolute error of labeling')
f.CLim=[0,250];
xlim([-4,6.5])
ylim([0,0.04])


%% Plot together
figure(5)
clf
set(gcf,'Position',[100,100,800,800]);

subplot(5,5,[2:5,7:10,12:15,17:20]);
f=densityScatterChart(sumPredictedLabels-trueSum,mae_FullLabels);
xlabel('Sum of predicted labels')
ylabel('Mean absolute error of labeling')
f.CLim=[0,200];

subplot(5,5,[22:25]);
histogram(sumPredictedLabels-trueSum,'NumBins',20);
set(gca,'TickDir','out','YDir','reverse')
axis off

subplot(5,5,[1:5:20]);
histogram(mae_FullLabels,'NumBins',20,'Orientation','horizontal');
set(gca,'XDir','reverse')
axis off