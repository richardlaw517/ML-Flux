%% Load Data
model = 'Simple'; % (Simple, UpperGly, Gly13C2H, GlyPPP, or CCM)

% Files can be reproduced via Data_Generation_and_Simulation
predFullFluxes = load(strcat('flux_pred_',model,'.dat'));
trueFreeFluxes = load(strcat('flux_test_',model,'.dat'));
kernelNet = load(strcat('KernelNet_',model,'.txt'));
kernelXch = load(strcat('KernelXch_',model,'.txt'));

%% Get full flux set for true free fluxes and normalize all fluxes
trueFullFluxes = zeros(size(predFullFluxes));

for i = 1:size(trueFreeFluxes,1)
    full_true_i = [kernelNet*trueFreeFluxes(i,1:size(kernelNet,2))';kernelXch*trueFreeFluxes(i,size(kernelNet,2)+1:end)'];
    trueFullFluxes(:,i) = (full_true_i./(full_true_i(1)))';
end; clear i full_true_i

predFullFluxes = predFullFluxes./predFullFluxes(:,1);

%% Mean errors

diffFluxes = predFullFluxes - trueFullFluxes; %

diff_avg = mean(diffFluxes,1,"omitmissing");
diff_std = std(diffFluxes,1,"omitmissing");

% GlyPPP
% freeNet = diff_avg([21,23:34]);
% fullNet = diff_avg([1:20,22,35:size(kernelNet,1)]);
% xch=diff_avg(size(kernelNet,1)+1:end);
% 
% GlyPPP without boundary transport fluxes
% freeNet = diff_avg([21,23]);
% fullNet = diff_avg([6:20,22]);
% xch=diff_avg(42:59);

% CCM
% freeNet = diff_avg([3,14:18,21,28,30:32,34,41,42,45,47:52]);
% fullNet = diff_avg([1:2,4:13,19:20,22:27,29,33,35:40,43,44,46]);
% xch=diff_avg(size(kernelNet,1)+1:end);

% CCM without boundary transport fluxes
% freeNet = diff_avg([21,28,30,34,41,42,45,48]);
% fullNet = diff_avg([4:13,19:20,22:27,29,33,35:40,43]);
% xch=diff_avg([56:65,71:82,85,87:95,97,100]);

Total=[freeNet,fullNet,xch];

Steps=[0,numel(freeNet),numel(freeNet)+numel(fullNet),numel(freeNet)+numel(fullNet)+numel(xch)];
[N,edges,bins] = histcounts(Total);  

figure(1);
clf
colors=lines;
hold on

for i=1:length(Steps)-1
    % histogram(Total(1+Steps(i):end),'BinEdges',[-Inf,-1.1:0.2:1.1,Inf],'FaceColor',colors(i,:));
    histogram(Total(1+Steps(i):end),'BinEdges',[-Inf,-0.11:0.02:0.11,Inf],'FaceColor',colors(i,:));
    hold on
end

% set(gca,'XLim',[-1.3,1.3])
set(gca,'XLim',[-.13,.13])
% set(gca,'YLim',[0,12])
legend(['Free net';'Full net';'Exchange'])

% xticks([-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1])
xticks([-0.11,-0.09,-0.07,-0.05,-0.03,-0.01,0.01,0.03,0.05,0.07,0.09,0.11])

set(gca,'TickDir','out');
