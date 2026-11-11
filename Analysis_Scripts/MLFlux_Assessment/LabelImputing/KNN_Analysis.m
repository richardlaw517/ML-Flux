%% Load Data  and pre-process (model dependent)
model = 'Simple'; % (Simple, UpperGly, Gly, GlyPPP, or CCM)
size = '100k'; % (100k or 1M)

% Source data files are either provided in KNN_Data or can also be reproduced via:
%   1) Generation of samples fluxes in Data_Generation_and_Simulation/SimulationScripts
%   2) Splitting into training, testing, and validation sets in
%      Training_Scripts/ at the random seed provided in the script
%   3  Prediction of labels by running the first block of Flux_Prediction/runML-Flux.py, which
%      returns the variable fullLabelSeet to be imported to predictedLabels here

label_test = readmatrix(strcat('KNN_Data/label_test_',model,'_',size,'.dat'));
label_dist = readmatrix(strcat('KNN_Data/label_dist_',model,'.txt'));
label_pred = readmatrix(strcat('KNN_Data/label_pred_',model,'.txt'));

labelsToKeep = readmatrix(strcat('LabelsToKeep/',model,'_Labels_to_Keep.xlsx'));
% NOTE: labelsToKeep is all the isotopomers physically possible, not
% about what's possible given the tracer (i.e.: B_M+2 and B_M+3 are still
% kept for the 1,0,0 A tracer in Simple despite it being impossible to
% label. This is because it's not obvious for every tracer
labelsToKeep = find(labelsToKeep == 1);

%% Calculate errors and distances
n = 10000; % 100 masks x 100 test labelings
maes = zeros(n,1);
maes_toKeep = zeros(n,1);
dist_tot = zeros(n,1);
for i = 1:n
    diff = mean(abs(label_test(mod(i-1,100)+1,:)-label_pred(i,:)));
    diff_toKeep = mean(abs(label_test(mod(i-1,100)+1,labelsToKeep)-label_pred(i,labelsToKeep)));
    
    maes(i) = diff;
    maes_toKeep(i) = diff_toKeep;
    dist_tot(i) = sum(label_dist(i,:));
end

%% Plot distance vs. MAE
x1 = dist_tot./3;
y1 = maes;
y1 = maes_toKeep;

figure('Renderer', 'painters', 'Position', [10 10 800 800])
clf

f=densityScatterChart(x1,y1);
f.CLim = [0,400];

% scatterhist(x1,y1,'Marker','.','NBins',[20,20])
% set(gca,'TickDir','out')
xlim([0,12E-3])
ylim([0,1.5E-4])

xlabel('Average Euclidean distance of 3 neighbors')
ylabel('Mean absolute error')