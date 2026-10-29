%% Preface: This code requires installation of the COBRA Toolbox
% doi.org/10.1038/s41596-018-0098-2

%% Load CCM_20240718  model and update constraints
load("CCM_Model.mat")

%% Sample fluxes
tic
options.nFiles=1;
options.nPointsPerFile=1000000;
[modelSampling,samples] = sampleCbModel(model,'SampledFluxes','ACHR',options);
toc

%% Conduct rejection sampling
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

%% Sampled fluxes
load("SampledNetFluxes_ACHR_20240718.dat") % Generated from above section

%% Exchange flux sampling
% NOTE: Since COBRA rejection sampling is focused purely on the
% distribution of net fluxes, exchange fluxes in the CCM model are sampled
% similary to other models (i.e., fluxinit, dependency-src)

xmlfile='CCM';
xlsname='CCM_Constraints_mat.xlsx';

% Create an EMU model
tic;
[model,memu]=modelinit(xmlfile,{'U'});
simulatefcn=str2func(xmlfile);
toc;

% Read constraints
free_net=zeros(size(model.kernel_net,2),1);
free_xch=zeros(size(model.kernel_xch,2),1);

[ineq,ineqStr]=xlsreadineq(xlsname,model,free_net,free_xch);
[eq,eqStr]=xlsreadeq(xlsname,model,free_net,free_xch);

% Randomize inital fluxes within constraints
n=17852; % Final size from ACHR sampling 
seed=0;
[free_net_set,free_xch_set]=fluxinit(model,ineq,eq,n,1,seed);

% Concanate to free_set (ML-Flux trained output)
free_set=[SampledNetFluxes_ACHR_20240718;free_xch_set]';