%% labeling simulation script

xmlfile='CCM';
xlsname='CCM_Constraints_1.xlsx';

%% create an EMU model
tic;
[model,memu]=modelinit(xmlfile,{'C'});
simulatefcn=str2func(xmlfile);
toc;

%% read constraints

free_net=zeros(size(model.kernel_net,2),1);
free_xch=zeros(size(model.kernel_xch,2),1);
[ineq,ineqStr]=xlsreadineq(xlsname,model,free_net,free_xch);
[eq,eqStr]=xlsreadeq(xlsname,model,free_net,free_xch);

%% randomize inital fluxes within constraints
total_fluxes = 1000000;
n = 1000;
seed = 1;
free_net_set = zeros(size(model.kernel_net,2),total_fluxes);
free_xch_set = zeros(size(model.kernel_xch,2),total_fluxes);
for iterPar = 1:total_fluxes/n
    tic
    [free_net_set(1:end,n*(iterPar-1)+1:n*iterPar),free_xch_set(1:end,n*(iterPar-1)+1:n*iterPar)]=par_fluxinit(model,ineq,eq,n,1,seed);
    disp(strcat(['Finished ',num2str(iterPar),' out of ',num2str(1000)]))
    toc
end
free_set=[free_net_set;free_xch_set]'; % OUTPUT for ML

 %% simulate isotope labeling 
tracers_glc = [1,1,1,1,1,1;
               1,0,0,0,0,0;
               0,1,0,0,0,0;
               1,1,0,0,0,0;
               0,0,0,0,0,1;
               1,0,0,0,0,1;
               0,0,1,0,0,0;
               0,0,0,0,1,0;
               0,0,0,1,0,0;
               0,0,0,0,1,1];

tracers_gln = [0,0,0,0,0;
               1,1,1,1,1];

for tracerConfig = 1:size(tracers_glc,1)*size(tracers_gln,1)                 % iterate through all tracers and combos
    % unlabeled tracers
    input.CO2_IN__C = isotopomervector([0],1);
    input.AC_IN__C = isotopomervector([0,0],1);
    input.OAA_IN__C = isotopomervector([0,0,0,0],1);
    input.Glu_IN__C = isotopomervector([0,0,0,0,0],1);

    % labeled tracers
    input.Gln_IN__C = isotopomervector(tracers_gln(mod(tracerConfig+1,2)+1,:),1); % alternate from unlabeled to labeled
    input.GLC_IN__C = isotopomervector(tracers_glc(round(tracerConfig/2),:),1); % 1,1,2,2,3,3,...9,9,10,10

    clear(char(simulatefcn))                            % reset the input metabolites, which is persistent in fcn

    for i=total_fluxes:-1:1                                        % run simulation with input metabolite
        simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
    end
    simulated_sets(tracerConfig)=mergestruct(simulated);          % INPUT for ML (one tracer)
end

% Extra simulation to replace the U-Glc U-Gln experiment with Unlabeled-Glc
% U-Gln in simulated_sets(2)
input.CO2_IN__C = isotopomervector([0],1);
input.AC_IN__C = isotopomervector([0,0],1);
input.OAA_IN__C = isotopomervector([0,0,0,0],1);
input.Glu_IN__C = isotopomervector([0,0,0,0,0],1);
input.Gln_IN__C = isotopomervector([1,1,1,1,1],1); % alternate from unlabeled to labeled
input.GLC_IN__C = isotopomervector([0,0,0,0,0,0],1); % 1,1,2,2,3,3,...9,9,10,10
clear(char(simulatefcn))                            % reset the input metabolites, which is persistent in fcn
for i=total_fluxes:-1:1                                        % run simulation with input metabolite
    simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
end
simulated_sets(2)=mergestruct(simulated);          % INPUT for ML (one tracer)

[~,simulated_set,lmid]=mergestruct2(simulated_sets);         % INPUT for ML (multiple tracers in struct)
singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)