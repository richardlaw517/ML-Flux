%% Load CCM_20240718  model and update constraints
load("CCM__Model.mat")

%% Sample fluxes
tic
options.nFiles=1;
options.nPointsPerFile=1000000;
[modelSampling,samples] = sampleCbModel_RL20240724(model,'SampledFluxes','ACHR',options);
toc

%% Conduct rejection sampling
load("SampledFluxes_1.mat") % Generated from above section
allACHRSamples = points;
reducedACHRSamples = allACHRSamples';

% Select a list of "key fluxes" to make the distribution more uniform
keyFluxes = {'ldh','pdh','ppck','me','mdh','akgdh'};

% Loop through each flux and cut off data to improve distribution
for keyFluxNum = 1:length(keyFluxes)
    
    dataToRemove = []; % Indices of data to remove
    fluxDistr = reducedACHRSamples(:,find(strcmp(model.rxns,keyFluxes{keyFluxNum}))); %#ok<FNDSB>
    
    num_bins = 20;
    max_count = 10000;
    
    [counts, edges] = histcounts(fluxDistr, num_bins);
   
    for i = 1:num_bins
        bin_indices = (fluxDistr >= edges(i)) & (fluxDistr < edges(i+1));
        
        num_points_in_bin = sum(bin_indices);
        
        if num_points_in_bin > max_count
            dataToRemove = [dataToRemove; randsample(find(bin_indices), num_points_in_bin-max_count)];
        end
    end

    % Remove points from current
    dataToRemove = unique(dataToRemove);
    reducedACHRSamples(dataToRemove,:) = [];
end