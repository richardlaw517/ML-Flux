%% Load ML-Flux Predicted fluxes
modelName = 'Simple'; % 'Simple','UpperGly','Glycolysis','GlyPPP',or 'CCM'

kernelNet = load(strcat('Kernels/KernelNet_',modelName,'.txt'));
kernelXch = load(strcat('Kernels/KernelXch_',modelName,'.txt'));
predictedFullFluxes = load(strcat('flux_pred_',modelName,'.dat'));
trueFreeFluxes = load(['flux_test_',modelName,'.dat']);

%% Get full flux set for true free fluxes and normalize all fluxes
trueFullFluxes = zeros(size(predictedFullFluxes));

for i = 1:size(trueFreeFluxes,1)
    full_true_i = [kernelNet*trueFreeFluxes(i,1:size(kernelNet,2))';kernelXch*trueFreeFluxes(i,size(kernelNet,2)+1:end)'];
    trueFullFluxes(i,:) = full_true_i./full_true_i(1);
end; clear i
predictedFullFluxes = predictedFullFluxes./predictedFullFluxes(:,1);

%% Get % errors of flux predictions
percentFluxError = (predictedFullFluxes - trueFullFluxes)./trueFullFluxes;

%% Find center percentiles and SD's for each flux
centralPercentile = 0.68269; % Corresponds to 1 stdev up and down

lowerPercentile = prctile(percentFluxError,50-centralPercentile*100/2);
upperPercentile = prctile(percentFluxError,50+centralPercentile*100/2);

normalStandardDev = (upperPercentile-lowerPercentile)'/2; 