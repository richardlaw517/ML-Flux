--
-%% Preface: This code requires installation of the COBRA Toolbox
% doi.org/10.1038/s41596-018-0098-2

%% Load CCM_20240718  model and update constraints
load("../Models and Constraints/mCM_model_constraints.mat")

%% reactions positions that will have knockouts
knockouts = [4:12,14:19,21:27,101:106]; 


%% loop through each reaction 
% set ub and lb to 0 and sample K fluxes with those constraints
K = 7000; 
% 7000 samples were generated for each of the 28 knockout reactions
% In order to maintain reasonable pyk>0 estimates for most knockouts:
% for reactions [8:14,16:28] 2000 samples were sampled from the constraints below and
% 5000 from the pyk>0 constraints commented out, for the remaining 8
% reactions all 7000 were from the constraints below.
 
allACHRSamples = [];

 for j= 1:28
  % for j= [8:14,16:28] % resampling for pyk>0 knockout constraints
    load('10012025_model.mat')
    model.lb(knockouts(1,1:4))=zeros(1,4); % upper glycolysis
    model.lb(8:11)=[-2; -2; -2; -2]; % lower glycolysis
    % Remove pyk<0 for pyk resampling
    %model.lb(8:10)=[-2; -2; -2; -2]; % lower glycolysis
    model.ub(12:18)=[1; 1; 1; 1; 1; 1; 1]; % PPP
    model.lb(13:17)=[-0.4; -0.4; -0.4; -0.4; -0.4]; %nonOxPPP
    model.ub(7)=2; %tpi
    model.lb(7)=-2; %tpi
    model.ub(101:105) = [1; 1; 1; 1; 1]; % ppc, me, oaadc, g6pdh, eda
    model.ub(92) = 2; % EX_DHAP_GLYC3P
    index= knockouts(1,j);
    model.lb(index)=0;
    model.ub(index)=0;
    
    % Identify blocked reactions
    blockedRxns=findBlockedReaction(model);
    [~, blockedIdx] = ismember(blockedRxns, model.rxns);

    % Sample fluxes
    tic
    options.nFiles=1;
    options.nPointsPerFile=K;
    [modelSampling,samples] = sampleCbModel(model,'SampledFluxes_KO_test','ACHR',options);
    toc
    % Load sampled fluxes
    load('SampledFluxes_KO_test_1.mat') %originally test_1

    % Add the zero rows back
    zero_row = zeros(1, K);
    for b = 1:length(blockedIdx)
    idx = blockedIdx(b);
    points = [points(1:idx-1,:); zero_row; points(idx:end,:)];
    end
    
    % Append to previous data
    allACHRSamples = [allACHRSamples, points];

 end
 



 
