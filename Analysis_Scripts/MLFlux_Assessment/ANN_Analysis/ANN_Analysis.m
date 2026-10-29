%% Load Data
model = 'CCM'; % Simple,UpperGly,Glycolysis,GlyPPP,CCM

% Source data files are either provided on the repo or can also be reproduced via:
%   1) Generation of samples fluxes in Data_Generation_and_Simulation/SimulationScripts
%   2) Splitting into training, testing, and validation sets in
%      Training_Scripts/ at the random seed provided in the script
%   3  Prediction by running through Flux_Prediction/runFluxPredictionFromFull.py, which
%      returns .dat (or .txt, interchangeably) files formatted for input

predFullFluxes = load(strcat('flux_pred_',model,'.dat'));
trueFreeFluxes = load(strcat('flux_test_',model,'.dat'));
kernelNet = load(strcat('KernelNet_',model,'.txt'));
kernelXch = load(strcat('KernelXch_',model,'.txt'));

%% Get full flux set for true free fluxes and normalize all fluxes
trueFullFluxes = zeros(size(predFullFluxes));

for i = 1:size(trueFreeFluxes,1)
    full_true_i = [kernelNet*trueFreeFluxes(i,1:size(kernelNet,2))';kernelXch*trueFreeFluxes(i,size(kernelNet,2)+1:end)'];
    trueFullFluxes(i,:) = (full_true_i./(full_true_i(1)))';
end; clear i full_true_i

predFullFluxes = predFullFluxes./predFullFluxes(:,1);

%% Mean errors (Fig.2f-g)
diffFluxes = predFullFluxes - trueFullFluxes; % ME

diff_avg = mean(diffFluxes,1,"omitmissing");
diff_std = std(diffFluxes,1,"omitmissing");

% GlyPPP
% freeNet = diff_avg([21,23:34]);
% fullNet = diff_avg([1:20,22,35:size(kernelNet,1)]);
% xch=diff_avg(size(kernelNet,1)+1:end);
% 
% GlyPPP without boundary transport fluxes
freeNet = diff_avg([21,23]);
fullNet = diff_avg([6:20,22]);
xch=diff_avg(42:59);

% CCM
% freeNet = diff_avg([3,14:18,21,28,30:32,34,41,42,45,47:52]);
% fullNet = diff_avg([1:2,4:13,19:20,22:27,29,33,35:40,43,44,46]);
% xch=diff_avg(size(kernelNet,1)+1:end);

% CCM without boundary transport fluxes
freeNet = diff_avg([21,28,30,34,41,42,45,48]);
fullNet = diff_avg([4:13,19:20,22:27,29,33,35:40,43]);
xch=diff_avg([56:65,71:82,85,87:95,97,100]);

Total=[freeNet,fullNet,xch];

Steps=[0,numel(freeNet),numel(freeNet)+numel(fullNet),numel(freeNet)+numel(fullNet)+numel(xch)];
[N,edges,bins] = histcounts(Total);  

% Plot stacked histogram
figure(1);
clf
colors=lines;
hold on

for i=1:length(Steps)-1
    histogram(Total(1+Steps(i):end),'BinEdges',[-Inf,-0.11:0.02:0.11,Inf],'FaceColor',colors(i,:));
    hold on
end

set(gca,'XLim',[-.13,.13])
legend(['Free net';'Full net';'Exchange'])

xticks([-0.11,-0.09,-0.07,-0.05,-0.03,-0.01,0.01,0.03,0.05,0.07,0.09,0.11])

set(gca,'TickDir','out');