%% Load data
% Source data files are either provided on the repo or can also be reproduced via:
%   1) Generation of samples fluxes in Data_Generation_and_Simulation/SimulationScripts
%   2) Splitting into training, testing, and validation sets in
%      Training_Scripts/ at the random seed provided in the script
%   3  Prediction by running through Flux_Prediction/runFluxPredictionFromFull.py, which
%      returns .dat (or .txt, interchangeably) files formatted for input
% Predicted fluxes for this are loaded as .mat files in the same folder

% Simple
load("dG_Simple.mat")

% Glycolysis
load("dG_Glycolysis.mat")

% CCM
load("dG_CCM.mat")

%% Process data (expand and normalize)
trueFullFluxes = zeros(size(trueFreeFluxes,1),size(predFullFluxes,2));

for i = 1:size(trueFreeFluxes,1)
    full_true_i = [kernelNet*trueFreeFluxes(i,1:size(kernelNet,2))';kernelXch*trueFreeFluxes(i,size(kernelNet,2)+1:end)'];
    trueFullFluxes(i,:) = (full_true_i./(full_true_i(1)))';
end; clear i full_true_i


predFullFluxes = predFullFluxes./predFullFluxes(:,1);
trueFullFluxes = repelem(trueFullFluxes,100,1);

% For Glycolysis only
% predFullFluxes = predFullFluxes(:,[1:20,22:end]);
% trueFullFluxes = trueFullFluxes(:,[1:20,22:end]);

% For CCM only (Includes 'unarrangeing' net fluxes)
rearrangeIdx = [1,2,3,21,22,23,24,25,26,27,28,29,5,11,12,13,14,15,30,52,51,31,32,33,34,35,36,37,38,39,16,17,40,18,41,42,43,44,45,46,47,48,19,20,49,7,8,50,4,6,9,10];
trueFullFluxes(:,1:52) = trueFullFluxes(:,rearrangeIdx);
predFullFluxes(:,1:52) = predFullFluxes(:,rearrangeIdx);


%% Calculate true and predicted dG
r = 8.314462618; % J/mol
t = 310; % K
numRxns = size(predFullFluxes,2)/2;

% For simple
% Pred
predNetFluxes = predFullFluxes(:,1:numRxns);
predXchFluxes = predFullFluxes(:,numRxns+1:end);
predForFluxes = abs(predXchFluxes + max(0,predNetFluxes));
predRevFluxes = abs(predXchFluxes - min(0,predNetFluxes));
dG_Pred = r .* t .* log(predRevFluxes./predForFluxes)./1000;
% True
trueNetFluxes = trueFullFluxes(:,1:numRxns);
trueXchFluxes = trueFullFluxes(:,numRxns+1:end);
trueForFluxes = trueXchFluxes + max(0,trueNetFluxes);
trueRevFluxes = trueXchFluxes - min(0,trueNetFluxes);
dG_True = r .* t .* log(trueRevFluxes./trueForFluxes)./1000;

% fluxesToEval = 2; % Simple
% fluxesToEval = 4:12; % Glycolysis
fluxesToEval = [5,7:11,13,23:29,37:38,40:43,45]; % CCM

%% For each flux, identify data in bins and make final dataset
% Skip this section for CCM as sampling is already distributed without need
% for binning

% binSize = 100;
% 
% dG_True_Binned = zeros(binSize*3,length(fluxesToEval));
% dG_Pred_Binned = zeros(binSize*3,length(fluxesToEval));
% small_dG_index = zeros(binSize,length(fluxesToEval));
% med_dG_index = zeros(binSize,length(fluxesToEval));
% large_dG_index = zeros(binSize,length(fluxesToEval));
% 
% for i = fluxesToEval
%     small_dG_index(:,i) = find((abs(dG_True(:,i))>=0)&(abs(dG_True(:,i))<0.5),binSize);
%     med_dG_index(:,i) = find((abs(dG_True(:,i))>=0.5)&(abs(dG_True(:,i))<1),binSize);
%     large_dG_index(:,i) = find((abs(dG_True(:,i))>=1)&(abs(dG_True(:,i))<5),binSize);
% 
%     dG_True_Binned(:,i) = dG_True([small_dG_index(:,i);med_dG_index(:,i);large_dG_index(:,i)],i);
%     dG_Pred_Binned(:,i) = dG_Pred([small_dG_index(:,i);med_dG_index(:,i);large_dG_index(:,i)],i);
% end; clear i;

%% Calculate the error
% deviation_dG = dG_Pred_Binned - dG_True_Binned; % For Simple and Glycolysis
deviation_dG = dG_Pred - dG_True; % For CCM

%% Get boxcharts of the apes

% reactionsList = {'v2'}; % Simple
% reactionsList = {'PGI','PFK','FBA','TPI','GAPD','PGK','PGM1','PGM2','ENO'}; % Glycolysis
reactionsList = {'pgi','fba','tpi','gapd','pgk','eno','ldh','gnd','rpi','rpe','tkt2','tkt1','tal','SBA','icdh','akgdh','sucd','fum','mdh','PYR_Ala','OGA_Glu'}; % CCM

xdata = repelem(reactionsList,1,size(deviation_dG,1))';
xdata = categorical(xdata);
xdata = reordercats(xdata,reactionsList);
ydata = deviation_dG(:,fluxesToEval);
ydata = ydata(:);

validIdx = ~isinf(ydata);
xdata = xdata(validIdx);
ydata = ydata(validIdx);

figure(1); clf; hold on
boxchartFig1 = boxchart(xdata,ydata);
set(boxchartFig1,'BoxWidth',0.8);
set(boxchartFig1,'MarkerStyle','none')
xlabel('Reaction')
ylabel('Deviation from true ∆G (predicted – true)')
ylim([-8,6])
ax = gca;
ax.YAxis.Scale ="linear";
set(gcf,'Position',[100,100,1000,500])
set(gca,'TickDir','out','FontSize',14,'LineWidth',1,'Xcolor','k','Ycolor','k')
% set(gca,'YTick',-0.08:0.02:0.06)