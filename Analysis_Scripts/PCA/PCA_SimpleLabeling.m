%% Load data
model = 'Simple';
fluxes = load(['fluxes_',model,'.txt']);
dataTable = readtable(['labeling_',model,'_PCA.csv']);
dataLabelsOnly = dataTable{:,2:end};

%% Process data for scaling, removal of features, etc.
means = mean(dataLabelsOnly,1);
sds = std(dataLabelsOnly,1);

% Uncentered, unscaled, removed M+0 and 0 variance values
labelsToKeep = sds~=0 & repmat([0,1,1,1],1,30); % Simple
featureNames = dataTable.Properties.VariableNames(2:end)';
featureNames = featureNames(labelsToKeep);
pca_Input = dataLabelsOnly(:,labelsToKeep);
%% Conduct PCA
[coeff,score,latent,tsquared,explained,mu] = pca(pca_Input);
covar = cov(pca_Input);
corr = corrcov(covar);

%% Find most important features in PCs 1 and 2
loadings = array2table(coeff(:,1:2));
loadings_abs = abs(loadings);
loadings = [featureNames,loadings];
weightedLoading = array2table(sum(abs(coeff).*explained'./100,2));
weightedLoading.Properties.VariableNames = "var3";

loadings_abs = [featureNames,loadings_abs,weightedLoading];
loadings.Properties.VariableNames = ["Isotopomer","PC1","PC2"];
loadings_abs.Properties.VariableNames = ["Isotopomer","PC1","PC2","Weighted loadings"];

%% Separate feature table by tracer experiment to find most important metabolites for each tracer
tracer = '';
metList = unique(cellfun(@(x) x(1:strfind(x,'M_')-2),featureNames,'UniformOutput',false));
loadingsOneTracer = loadings_abs(contains(table2cell(loadings_abs(:,1)),tracer),:);
loadingsOneTracerSummedMets = zeros(length(metList),size(loadingsOneTracer,2)-1);

for i = 1:length(metList)
    loadingsOneTracerSummedMets(i,:) = table2array(sum(loadingsOneTracer(contains(table2cell(loadingsOneTracer(:,1)),metList(i)),2:end),1));
end
loadingsOneTracerSummedMets = array2table(loadingsOneTracerSummedMets);
loadingsOneTracerSummedMets = [metList,loadingsOneTracerSummedMets];
loadingsOneTracerSummedMets.Properties.VariableNames = ["Isotopomer","PC1","PC2","Weighted loadings"];

%% Scree plot
screeData = [explained(1:2);sum(explained(3:end))];
cummulativeExplanation = cumsum(screeData);

figure(1); clf; hold on;
plot_Scree = bar(screeData);
screeSum = plot(cummulativeExplanation,'LineWidth',1);
xlim([0.5,length(screeData)+0.5])
ylim([0,100])
xticks(1:length(screeData))
xticklabels({'1','2','All other PCs'})
xlabel('Principal component')
ylabel('Percent of variance captured by component')
set(gca,'TickDir','out')
legend({'Variance captured in PC','Cummulative variance captured'},'Location','east')
set(gca,'TickDir','out','FontSize',14,'LineWidth',1,'Xcolor','k','Ycolor','k')

%% 2D PC plot
pc_x = 1;
pc_y = 2;
fluxNum = 3;

data_PC_x = pca_Input*coeff(:,pc_x);
data_PC2_y = pca_Input*coeff(:,pc_y);
data_flux = fluxes(:,fluxNum);

figure(1); clf; hold on
plot2_PCA = scatter(data_PC_x,data_PC2_y,[],data_flux,'.');

xlabel(['PC',num2str(pc_x),' (',num2str(explained(pc_x)),'%)'])
ylabel(['PC',num2str(pc_y),' (',num2str(explained(pc_y)),'%)'])
title('v_2_r as a function of PCs')
set(gca,'TickDir','out')

set(gca,'ColorScale','log')
clim([1,10])
colorbar('Ticks', 1:1:10,'TickLabels', ...
    {'<=1','2','3','4','5','6','7','8','9','>=10'},'TickDirection','out');

%% PC loading bar plot of weighted absolute loadings sorted
figure(2); clf; hold on
set(gcf,'Position',[200,200,1500,500])
[sortedWLoadings, sortOrder] = sort(table2array(weightedLoading),'descend');
sortedFeatures = featureNames(sortOrder);
barLoading_weighted = bar(sortedFeatures,sortedWLoadings);
set(gca,'TickLabelInterpreter','none','TickDir','out')
xlabel('Isotopomer (not including M+0s or any 0 variance isotopomers)')
ylabel('PC importance weighted loadings of each isotopomer')
set(gca,'TickDir','out','FontSize',14,'LineWidth',1,'Xcolor','k','Ycolor','k')
