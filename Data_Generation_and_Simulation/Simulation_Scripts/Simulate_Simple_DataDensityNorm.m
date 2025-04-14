%% Setup: Read/generate model and constraints
% Directory to model and constraints file
xmlfile='Simple';
xlsname='Simple_Constraints.xlsx';

% Generation of EMU model
[model,~]=modelinit(xmlfile,{'U'});
simulatefcn=str2func(xmlfile);

% Generation of constrains variables
free_net=zeros(size(model.kernel_net,2),1);
free_xch=zeros(size(model.kernel_xch,2),1);
[ineq,~]=xlsreadineq(xlsname,model,free_net,free_xch);
[eq,~]=xlsreadeq(xlsname,model,free_net,free_xch);

clear xmlfile xlsname

%% Set parameters for flux sampling
seed = 0; % Starting random number seed
numFluxes = 5000; % Number of flux-label datapoints to generate
nFirstBatch = 100; % Size of first batch
sampleBatchSize = 10000; % Size of fluxes to sample from fluxinit
k = 10; % Number of neighbors

%% Define tracers to simulate
allTracers  = [1,0,0;
               0,1,0;
               1,1,0;
               0,0,1;
               1,0,1;
               0,1,1];

%% Generate large batch of flux-label points to sample from
[free_set_batch,labeling_batch] = generateFluxLabelBatch(model,ineq,eq,sampleBatchSize,seed,simulatefcn,allTracers);

%% Generation of first n flux-label points and initiation of final matrices
% Sample first flux-label point
free_set_i = free_set_batch(1:nFirstBatch,:);
labeling_i = labeling_batch(1:nFirstBatch,:);

% Create final flux and labeling matrices and add first nFirstBatch points
free_set = zeros(numFluxes,size(free_set_i,2)); % OUTPUT of ML-Flux
free_set(1:nFirstBatch,:) = free_set_i;
singlemat = zeros(numFluxes,size(labeling_i,2));
singlemat(1:nFirstBatch,:) = labeling_i; % INPUT of ML-Flux

% Make a matrix to save all final distances and probabilities
allMeanMinKDists = zeros(numFluxes,1);

% Calculate and add first nFirstBatch distances
dist = pdist2(singlemat(1:nFirstBatch,:),singlemat(1:nFirstBatch,:));
dist(dist==0) = NaN;
allMeanMinKDists(1:nFirstBatch) = mean(mink(dist,k));

%% Generation of remaining numFluxes points with DDN
sampleBatchSampleNum = nFirstBatch;
for fluxNum_i = nFirstBatch+1:numFluxes % Begin after nFirstBatch
    sampleBatchSampleNum =  sampleBatchSampleNum+1;
    if sampleBatchSampleNum > sampleBatchSize
        seed = seed+1;
        [free_set_batch,labeling_batch] = generateFluxLabelBatch(model,ineq,eq,sampleBatchSize,seed,simulatefcn,allTracers);
        sampleBatchSampleNum = 1;
    end
    
    % Sample a flux-label point
    free_set_i = free_set_batch(sampleBatchSampleNum,:);
    labeling_i = labeling_batch(sampleBatchSampleNum,:);

    % Calculate distance of points
    dist = pdist2(singlemat(1:fluxNum_i-1,:),labeling_i);
    meanMinKDist = mean(mink(dist,k));

    % Get a distribution of current accepted distances
    prob = size(allMeanMinKDists(allMeanMinKDists(1:fluxNum_i)<meanMinKDist),1)/fluxNum_i;
    randomR = rand;

    while randomR >= prob
        sampleBatchSampleNum =  sampleBatchSampleNum+1;
        if sampleBatchSampleNum > sampleBatchSize
            seed = seed+1;
            [free_set_batch,labeling_batch] = generateFluxLabelBatch(model,ineq,eq,sampleBatchSize,seed,simulatefcn,allTracers);
            sampleBatchSampleNum = 1;
        end
        free_set_i = free_set_batch(sampleBatchSampleNum,:);
        labeling_i = labeling_batch(sampleBatchSampleNum,:);
        dist = pdist2(singlemat(1:fluxNum_i-1,:),labeling_i);
        meanMinKDist = mean(mink(dist,k));
        prob = size(allMeanMinKDists(allMeanMinKDists(1:fluxNum_i)<meanMinKDist),1)/fluxNum_i;
        prob2 = size(allMeanMinKDists(allMeanMinKDists(1:fluxNum_i).^0.5<meanMinKDist.^0.5),1)/(fluxNum_i-1);
        randomR = rand;
    end
    free_set(fluxNum_i,:) = free_set_i;
    singlemat(fluxNum_i,:) = labeling_i;
    allMeanMinKDists(fluxNum_i) = meanMinKDist;
end

%% Functions
function [free_set_batch,labeling_batch] = generateFluxLabelBatch(model,ineq,eq,sampleBatchSize,seed,simulatefcn,allTracers)
    [free_net_set,free_xch_set]=fluxinit_logUniXch(model,ineq,eq,sampleBatchSize,1,seed);
    free_set_batch=[free_net_set;free_xch_set]'; % OUTPUT for ML
    
    for tracerNum = 1:size(allTracers,1)
        input.molA__U = isotopomervector(allTracers(tracerNum,:),1);
        clear(char(simulatefcn)) % Reset the input metabolite, which is persistent in fcn
        for i=1:sampleBatchSize
            simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
        end
        simulated_sets(tracerNum) = mergestruct(simulated); 
    end
    [~,simulated_set,~] = mergestruct2(simulated_sets);
    labeling_batch = struct2mat(simulated_set); % matrix of labeling

end