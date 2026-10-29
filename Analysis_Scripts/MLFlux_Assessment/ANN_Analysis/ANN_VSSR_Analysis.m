%% Load ML-Flux Predicted fluxes
modelName = 'CCM'; % 'Simple','UpperGly','Glycolysis','GlyPPP',or 'CCM'

% Source data files are provided or can also be reproduced via Data_Generation_and_Simulation
predictedFullFluxes = load(strcat('flux_pred_',modelName,'.dat'));
kernelNet = load(strcat('KernelNet_',modelName,'.txt'));
kernelXch = load(strcat('KernelXch_',modelName,'.txt'));
trueFreeFluxes = load(strcat('flux_test_',modelName,'.dat'));

% True labeling and fluxes
trueFullLabeling = load(['label_test_',modelName,'.dat']);

%% Get full flux set for true free fluxes and normalize all fluxes
trueFullFluxes = zeros(size(predictedFullFluxes));

for i = 1:size(trueFreeFluxes,1)
    full_true_i = [kernelNet*trueFreeFluxes(i,1:size(kernelNet,2))';kernelXch*trueFreeFluxes(i,size(kernelNet,2)+1:end)'];
    trueFullFluxes(i,:) = full_true_i./full_true_i(1);
end; clear i
predictedFullFluxes = predictedFullFluxes./predictedFullFluxes(:,1);

%% Get the free_set of fluxes (input for labeling simulation)
% NOTE: Only run if simulating labeling (not necessary for CCM, see below)
% if strcmp(modelName,'Simple')
%     free_net_index = 6:7; free_xch_index = 9:12;
% elseif strcmp(modelName,'UpperGly')
%     free_net_index = 14:21; free_xch_index = 25:34;
% elseif strcmp(modelName,'Glycolysis')
%     free_net_index = 14:21; free_xch_index = 25:34;
% elseif strcmp(modelName,'GlyPPP')
%     free_net_index = [21,23:34]; free_xch_index = [42:59,71:72];
% elseif strcmp(modelName,'CCM')
%     free_net_index = [3,14:18,21,28,30:32,34,41:42,45,47:52]; free_xch_index = [56:65,71:82,85,87:95,97,100];
% end
% 
% free_net_set = predictedFullFluxes(:,free_net_index)';
% free_xch_set = predictedFullFluxes(:,free_xch_index)';
% free_set=[free_net_set;free_xch_set]';
% 
% free_net_set_True = trueFullFluxes(:,free_net_index)';
% free_xch_set_True = trueFullFluxes(:,free_xch_index)';
% free_set_True=[free_net_set;free_xch_set]';

%% Create EMU model
% NOTE: Will require access to 'src' in Data_Generation_and_Simulation
% if strcmp(modelName,'Simple')
%     xmlfile='Simple';
% elseif strcmp(modelName,'UpperGly')
%     xmlfile='Glycolysis';
% elseif strcmp(modelName,'Glycolysis')
%     xmlfile='Glycolysis';
% elseif strcmp(modelName,'GlyPPP')
%     xmlfile='GlyPPP';
% elseif strcmp(modelName,'CCM')
%     xmlfile='CCM';
% end
% 
% [model,memu]=modelinit(xmlfile,{'U'});
% simulatefcn=str2func(xmlfile);

%% Simulate isotope labeling
% NOTE: This may take a large amount of time depending on computing power.
% Labeling simulations for the main (CCM) model are provided and can be
% loaded in the next block
% if strcmp(modelName,'Simple')
%     m1=3;                                                   % # of atoms in input metabolite 1; may need multiple inpute mets
%     for j1=2^m1-2:-1:1                                      % generate all nontrivial tracers
%         alltracers(j1,:)=de2bi(j1,m1);
%         input.molA__U=isotopomervector(alltracers(j1,:),1);
%         clear(char(simulatefcn))% reset the input metabolite, which is persistent in fcn
%         for i=1:size(free_set,1)
%             simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
%         end
%         simulated_sets(j1)=mergestruct(simulated);
%     end
%     [~,simulated_set,lmid]=mergestruct2(simulated_sets);
%     singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)
% end
% if strcmp(modelName,'UpperGly')
%     tracers = [1,1,0,0,0,0,0,0,0,0,1,0,0];                  % List out tracer configs used
%     for j1=1:size(tracers,1)                                % iterate through all tracers
%         input.GLC_IN__C=isotopomervector(tracers(j1,1:6),1);
%         input.GLC_IN__H=isotopomervector(tracers(j1,:),1);
%         input.H_IN__H=isotopomervector([0],1);
%         input.NADH_IN__H=isotopomervector([0],1);
%         clear(char(simulatefcn))% reset the input metabolites, which is persistent in fcn
%         for i=1:size(free_set,1)
%             simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
%         end
%         simulated_sets(j1)=mergestruct(simulated);
%     end
%     [~,simulated_set,lmid]=mergestruct2(simulated_sets);
%     singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)
% end
% if strcmp(modelName,'Glycolysis')
%     tracers = [1,1,0,0,0,0,0,0,0,0,1,0,0];                  % List out tracer configs used
%     for j1=1:size(tracers,1)                                % iterate through all tracers
%         input.GLC_IN__C=isotopomervector(tracers(j1,1:6),1);
%         input.GLC_IN__H=isotopomervector(tracers(j1,:),1);
%         input.H_IN__H=isotopomervector([0],1);
%         input.NADH_IN__H=isotopomervector([0],1);
%         clear(char(simulatefcn))% reset the input metabolites, which is persistent in fcn
%         for i=1:size(free_set,1)
%             simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
%         end
%         simulated_sets(j1)=mergestruct(simulated);
%     end
%     [~,simulated_set,lmid]=mergestruct2(simulated_sets);
%     singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)
% end
% if strcmp(modelName,'GlyPPP')
%     input.CO2_IN__U=isotopomervector([0],1);
% 
%     tracers=[1,1,1,0,0,0,0,0,0,0,0,0,0;
%              1,1,0,0,0,0,0,0,0,0,0,0,0;
%              1,0,0,0,0,0,0,0,0,0,0,0,0;
%              0,1,1,1,1,1,0,0,0,0,0,0,0;
%              0,1,1,0,0,0,0,0,0,0,0,0,0;
%              0,1,0,0,0,0,0,0,0,0,0,0,0;
%              0,0,1,1,0,0,0,0,0,0,0,0,0;
%              0,0,1,0,0,0,0,0,0,0,0,0,0;
%              0,0,0,1,1,1,0,0,0,0,0,0,0;
%              0,0,0,1,1,0,0,0,0,0,0,0,0;
%              0,0,0,1,0,0,0,0,0,0,0,0,0;
%              0,0,0,0,1,0,0,0,0,0,0,0,0;
%              0,0,0,0,0,1,0,0,0,0,0,0,0;
%              1,1,1,1,1,1,0,0,0,0,0,0,0;
%              0,1,0,0,1,1,0,0,0,0,0,0,0;
%              1,0,0,0,1,1,0,0,0,0,0,0,0;
%              1,1,0,1,0,0,0,0,0,0,0,0,0;
%              0,1,0,1,0,0,0,0,0,0,0,0,0;
%              1,0,0,1,0,0,0,0,0,0,0,0,0;
%              0,0,0,1,0,1,0,0,0,0,0,0,0;
%              0,1,0,0,1,0,0,0,0,0,0,0,0;
%              1,0,1,0,0,0,0,0,0,0,0,0,0;
%              0,0,0,0,1,1,0,0,0,0,0,0,0;
%              1,0,0,0,0,1,0,0,0,0,0,0,0;];%List out tracer configs used
% 
%     for tracerConfig = 1:size(tracers,1)
%         tic
%         input.GLC_IN__U=isotopomervector(tracers(tracerConfig,:),1);
%         clear(char(simulatefcn))    % reset the input metabolite, which is persistent in fcn
% 
%         parfor i=1:size(free_set,1)                                        % run simulation with input metabolite
%             simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
%             simulated_True(i)=simulatefcn(free_net_set_True(:,i),free_xch_set_True(:,i),input);
%         end
% 
%         simulated_sets(tracerConfig)=mergestruct(simulated);          % INPUT for ML (one tracer)
%         simulated_sets_True(tracerConfig)=mergestruct(simulated_True);
%         disp(strcat(['Finished ',num2str(tracerConfig),' fluxes out of ',num2str(size(tracers,1))]))
%         toc
%     end
%     [~,simulated_set,lmid]=mergestruct2(simulated_sets); 
%     singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)
%     [~,simulated_set_True,lmid_True]=mergestruct2(simulated_sets_True); 
%     trueFullLabeling=struct2mat(simulated_set_True);                    % INPUT for ML (multiple tracers in a single matrix)
% 
% end
% if strcmp(modelName,'CCM')
%     tracers_glc = [1,1,1,1,1,1;
%                    1,0,0,0,0,0;
%                    0,1,0,0,0,0;
%                    1,1,0,0,0,0;
%                    0,0,0,0,0,1;
%                    1,0,0,0,0,1;
%                    0,0,1,0,0,0;
%                    0,0,0,0,1,0;
%                    0,0,0,1,0,0;
%                    0,0,0,0,1,1];
% 
%     tracers_gln = [0,0,0,0,0;
%                    1,1,1,1,1];
% 
%     for tracerConfig = 1:size(tracers_glc,1)*size(tracers_gln,1)                 % iterate through all tracers and combos
%         tic
%         % unlabeled tracers
%         input.CO2_IN__C = isotopomervector([0],1);
%         input.AC_IN__C = isotopomervector([0,0],1);
%         input.OAA_IN__C = isotopomervector([0,0,0,0],1);
%         input.Glu_IN__C = isotopomervector([0,0,0,0,0],1);
% 
%         % labeled tracers
%         input.Gln_IN__C = isotopomervector(tracers_gln(mod(tracerConfig+1,2)+1,:),1); % alternate from unlabeled to labeled
%         input.GLC_IN__C = isotopomervector(tracers_glc(round(tracerConfig/2),:),1); % 1,1,2,2,3,3,...9,9,10,10
% 
%         clear(char(simulatefcn))                            % reset the input metabolites, which is persistent in fcn
% 
%         parfor i=1:size(free_set,1)                                        % run simulation with input metabolite
%             simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
%             simulated_True(i)=simulatefcn(free_net_set_True(:,i),free_xch_set_True(:,i),input);
%         end
%         simulated_sets(tracerConfig)=mergestruct(simulated);
%         simulated_sets_True(tracerConfig)=mergestruct(simulated_True);
%         disp(strcat(['Finished ',num2str(tracerConfig),' fluxes out of ',num2str(20)]))
%         toc
%     end
% 
%     % Extra simulation to replace the U-Glc U-Gln experiment with Unlabeled-Glc
%     % U-Gln in simulated_sets(2)
%     input.CO2_IN__C = isotopomervector([0],1);
%     input.AC_IN__C = isotopomervector([0,0],1);
%     input.OAA_IN__C = isotopomervector([0,0,0,0],1);
%     input.Glu_IN__C = isotopomervector([0,0,0,0,0],1);
%     input.Gln_IN__C = isotopomervector([1,1,1,1,1],1); % alternate from unlabeled to labeled
%     input.GLC_IN__C = isotopomervector([0,0,0,0,0,0],1); % 1,1,2,2,3,3,...9,9,10,10
%     clear(char(simulatefcn))                            % reset the input metabolites, which is persistent in fcn
%     parfor i=1:size(free_set,1)                                       % run simulation with input metabolite
%         simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
%         simulated_True(i)=simulatefcn(free_net_set_True(:,i),free_xch_set_True(:,i),input);
%     end
%     simulated_sets(2)=mergestruct(simulated);
%     simulated_sets_True(2)=mergestruct(simulated_True);
% 
%     [~,simulated_set,lmid]=mergestruct2(simulated_sets);         % INPUT for ML (multiple tracers in struct)
%     singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)
%     [~,simulated_set_True,lmid_True]=mergestruct2(simulated_sets_True); 
%     trueFullLabeling=struct2mat(simulated_set_True);         % INPUT for ML (multiple tracers in a single matrix)
% end
% 
% singlemat = full(singlemat);


%% Load labeling data (for CCM simulations from mfapy)
singlemat = load("CCM_LabelingSimulatedFromPredictedFluxes.dat")';

%% Calculate VSSR for predicted simulated labeling versus true labeling and MAPE fluxes
ssr_pred_full = sum(((singlemat - trueFullLabeling).^2),2);
diffLabels = singlemat - trueFullLabeling;
abs_diffLabels = abs(singlemat - trueFullLabeling);

% Remove transport fluxes for GlyPPP or CCM models
% removeList = [1:5,24:41,51,60:size(trueFullFluxes,2)]; % GlyPPP
removeList = [1:4,6:18,20,21,53:55,66:70,83:84,86,88,96,98:99,101:104]; % CCM
removed_trueFullFluxes = trueFullFluxes;
removed_predictedFullFluxes = predictedFullFluxes;
removed_trueFullFluxes(:,removeList) = [];
removed_predictedFullFluxes(:,removeList) = [];

mape_allFluxes = mean(abs(round(removed_trueFullFluxes-removed_predictedFullFluxes,5)./removed_trueFullFluxes),2,'omitnan');

%% Get only net or exchange fluxes
% GlyPPP
net_index = 1:18;
xch_index = 19:35;

% CCM
net_index = 1:33;
xch_index = [35,37:41,43,48:54,58,59,61:65,];

mae_allFluxes = mean(abs(round(removed_trueFullFluxes(:,net_index)-removed_predictedFullFluxes(:,net_index),5)),2,'omitnan');
mape_allFluxes = mean(abs(round(removed_trueFullFluxes(:,net_index)-removed_predictedFullFluxes(:,net_index),5)./removed_trueFullFluxes(:,net_index)),2,'omitnan');
mape_vals = abs(round(removed_trueFullFluxes(:,net_index) - removed_predictedFullFluxes(:,net_index), 5) ./ removed_trueFullFluxes(:,net_index));
mape_vals(isinf(mape_vals)) = NaN; % Replace inf values with NaN
mape_allFluxes = mean(mape_vals, 2, 'omitnan');

%% Histograms for MAPE and SSR (SuppFig.2)
figure(1); clf; hold on
tiledlayout(2,1); nexttile;
histogram(mape_allFluxes*100,logspace(0,4,20))
set(gca,'XScale','log','Xlim',[1E-1,1E3],'TickDir','out')
title('MAPE')
nexttile;
histogram(ssr_pred_full,logspace(-1,3,20))
set(gca,'XScale','log','Xlim',[1E-2,1E3],'TickDir','out')
title('SSR')

%% Density in each flux prediction (SuppFig.2)
figure(2); clf;
figure2 = densityScatterChart(mape_allFluxes*100,ssr_pred_full);
xlim([1E0,1E4])
ylim([1E-1,1E3])
figure2.CLim = [0,1000];
figure2.ColorScale = 'linear';
set(gca,'yscale','log','xscale','log')
set(gca,'TickDir','out','FontSize',14,'LineWidth',1)
xlabel('Mean absolute percent error of fluxes')
ylabel('SSR of labeling from predicted fluxes')
set(gcf,'Position',[100,100,550,500])
