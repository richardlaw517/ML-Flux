%% Load Data
model = "Simple"; % Change this to Simple, UpperGly, Glycolysis, GlyPPP, or CCM

predictedFullFluxes = load(strcat('PredictedData/',model,'_FullPredictedFluxes_Inpaint_Unique.txt'));
predictedFullFluxesFromTrueLabels = load(strcat('PredictedData/',model,'_FullPredictedFluxesFromTrueLabels.txt'));
trueFreeFluxes = load(strcat('TestingData/',model,'_Fluxes_100_Seed0.txt'));
kernelNet = load(strcat('TestingModels/FluxPrediction/Kernels/KernelNet_',model,'.txt'));
kernelXch = load(strcat('TestingModels/FluxPrediction/Kernels/KernelXch_',model,'.txt'));

%% Get full flux set for true free fluxes and normalize all fluxes
trueFullFluxes = zeros(size(predictedFullFluxes));
trueFullFluxesFromTrueLabels = zeros(size(trueFreeFluxes,1),size(trueFullFluxes,2));

for i = 1:size(trueFreeFluxes,1)
    full_true_i = [kernelNet*trueFreeFluxes(i,1:size(kernelNet,2))';kernelXch*trueFreeFluxes(i,size(kernelNet,2)+1:end)'];
    trueFullFluxes(100*(i-1)+1:100*(i),:) = repmat((full_true_i./full_true_i(1))',100,1);
    % trueFullFluxes(31*(i-1)+1:31*(i),:) = repmat((full_true_i)',31,1);
    trueFullFluxesFromTrueLabels(i,:) = full_true_i;
end; clear i full_true_i

predictedFullFluxes = predictedFullFluxes./predictedFullFluxes(:,1);
predictedFullFluxesFromTrueLabels = predictedFullFluxesFromTrueLabels./predictedFullFluxesFromTrueLabels(:,1);

%% Get mean differences between predicted true and predicted scratched
diff_TrueScratched = predictedFullFluxes - repelem(predictedFullFluxesFromTrueLabels,100,1);
% diff_TrueScratched = predictedFullFluxes - repelem(predictedFullFluxesFromTrueLabels,31,1);
diff_avg_TrueScratched = mean(diff_TrueScratched,1);
diff_std_TrueScratched = std(diff_TrueScratched,1);

%% Get errors of predicted fluxes
% Absolute error
diffFluxesFromScratched_ae = abs(predictedFullFluxes - trueFullFluxes);
diffFluxesFromTrueLabels_ae = abs(predictedFullFluxesFromTrueLabels - trueFullFluxesFromTrueLabels);

% Absolute percent error
diffFluxesFromScratched_ape = abs((predictedFullFluxes - trueFullFluxes)./trueFullFluxes);
diffFluxesFromTrueLabels_ape = abs((predictedFullFluxesFromTrueLabels - trueFullFluxesFromTrueLabels)./trueFullFluxesFromTrueLabels);

diff_avg_FromScratched_ae = mean(diffFluxesFromScratched_ae,1,"omitmissing");
diff_std_FromScratched_ae = std(diffFluxesFromScratched_ae,1,"omitmissing");

diff_avg_FromTrue_ae = mean(diffFluxesFromTrueLabels_ae,1,"omitmissing");
diff_std_FromTrue_ae = std(diffFluxesFromTrueLabels_ae,1,"omitmissing");

diff_avg_FromScratched_ape = mean(diffFluxesFromScratched_ape,1,"omitmissing");
diff_std_FromScratched_ape = std(diffFluxesFromScratched_ape,1,"omitmissing");

diff_avg_FromTrue_ape = mean(diffFluxesFromTrueLabels_ape,1,"omitmissing");
diff_std_FromTrue_ape = std(diffFluxesFromTrueLabels_ape,1,"omitmissing");

diffFluxesFromScratched = diffFluxesFromScratched_ape;
diffFluxesFromTrueLabels = diffFluxesFromTrueLabels_ape;

diff_avg_FromScratched = diff_avg_FromScratched_ae;
diff_std_FromScratched = diff_std_FromScratched_ae;

diff_avg_FromTrue = diff_avg_FromTrue_ae;
diff_std_FromTrue = diff_std_FromTrue_ae;

%% Get boxcharts of the errors

reactionsList = load("FullReactionNames.mat",strcat(model));
reactionsList = reactionsList.(strcat(model));
if model == "Simple"
    removeList = [1,8,10:14]; % Removes the basis flux, transport fluxes, and irrelevant exchange fluxes
elseif model == "UpperGly"
    removeList = 1:3;
elseif model == "Glycolysis"
    removeList = [1:24,34:41];
elseif model == "GlyPPP"
    removeList = [1:5,24:41,43,51,59:size(reactionsList,2)];
elseif model == "CCM"
    removeList = [1:4,5:18,20,21,53:55,66:70,83:84,86,88,96,98:99,101:104];
end
reactionsList(removeList) = [];

xdata = [repelem(reactionsList,1,100)';repelem(reactionsList,1,10000)'];
xdata = categorical(xdata);
xdata = reordercats(xdata,reactionsList);
diffFluxesFromTrueLabels_Box = diffFluxesFromTrueLabels;
diffFluxesFromTrueLabels_Box(:,removeList) = [];
diffFluxesFromScratched_Box = diffFluxesFromScratched;
diffFluxesFromScratched_Box(:,removeList) = [];
ydata = [diffFluxesFromTrueLabels_Box(:);diffFluxesFromScratched_Box(:)];
ydata(ydata<1E-7) = NaN;
type = [repelem({'All labels'},numel(diffFluxesFromTrueLabels_Box))';repelem({'Scratched labels'},numel(diffFluxesFromScratched_Box))'];

figure(2)
clf
hold on

boxchartFig = boxchart(xdata,ydata,'GroupByColor',type);
set(boxchartFig,'BoxWidth',0.8);
set(boxchartFig,'MarkerStyle','none')
legend
ylim([1E-7,1E3])
ax = gca;
ax.YAxis.Scale ="log";
set(gcf,'Position',[100,100,1200,600])
set(gca,'TickDir','out')

%% Get boxcharts of the errors (Split by magnitude - cutoff 0.05)
avg_true = mean(abs(trueFullFluxesFromTrueLabels));
mape_list = find(avg_true>=0.05);
mae_list = find(avg_true<0.05);

reactionsList = load("FullReactionNames.mat",strcat(model));
reactionsList = reactionsList.(strcat(model));
if model == "Simple"
    removeList = [1,8,10:14]; % Removes the basis flux, transport fluxes, and irrelevant exchange fluxes
elseif model == "UpperGly"
    removeList = 1:3;
elseif model == "Glycolysis"
    removeList = [1:24,34:41];
elseif model == "GlyPPP"
    removeList = [1:5,24:41,43,51,59:size(reactionsList,2)];
elseif model == "CCM"
    removeList = [1:4,6:18,20,21,53:55,66:70,83:84,86,88,96,98:99,101:104];
end
reactionsList = reactionsList(setdiff(mae_list,removeList));

xdata = [repelem(reactionsList,1,100)';repelem(reactionsList,1,10000)'];
xdata = categorical(xdata);
xdata = reordercats(xdata,reactionsList);
diffFluxesFromTrueLabels_Box = diffFluxesFromTrueLabels_ae(:,setdiff(mae_list,removeList));
diffFluxesFromScratched_Box = diffFluxesFromScratched_ae(:,setdiff(mae_list,removeList));
ydata = [diffFluxesFromTrueLabels_Box(:);diffFluxesFromScratched_Box(:)];
ydata(ydata<1E-7) = NaN;
type = [repelem({'All labels'},numel(diffFluxesFromTrueLabels_Box))';repelem({'Scratched labels'},numel(diffFluxesFromScratched_Box))'];

figure(2)
clf
hold on

boxchartFig = boxchart(xdata,ydata,'GroupByColor',type);
set(boxchartFig,'BoxWidth',0.8);
set(boxchartFig,'MarkerStyle','none')
legend
ylim([1E-7,1E3])
% ylim([0,0.6])
ax = gca;
ax.YAxis.Scale ="log";
set(gcf,'Position',[100,100,1200,600])
set(gca,'TickDir','out')


%% Stacked histogram from full and scratched (For GlyPPP and CCM model)

% All non-GlyPPP non_CCM models
% A=diff_avg_FromTrue;
% B=diff_avg_FromScratched;

% Only the non-EX fluxes GlyPPP
% A=[diff_avg_FromTrue(6:23),diff_avg_FromTrue(42:59),];
% B=[diff_avg_FromScratched(6:23),diff_avg_FromScratched(42:59),];

% Only the non-EX fluxes CCM
A=diff_avg_FromTrue([5,19,22:52,56:65,71:82,85,87:95,97,100]);
B=diff_avg_FromScratched([5,19,22:52,56:65,71:82,85,87:95,97,100]);


% Plot stacked histogram   
figure(3);
clf
hold on

histogram(A,'BinEdges',[0,0:0.05:1.0,inf],'facealpha',.7,'edgecolor','none')
histogram(B,'BinEdges',[0,0:0.05:1.0,inf],'facealpha',.3,'edgecolor','none')
xlabel('Mean absolute error')
ylabel('Number of fluxes')
set(gca,'TickDir','out')