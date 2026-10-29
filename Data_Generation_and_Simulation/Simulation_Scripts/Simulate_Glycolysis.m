%% labeling simulation script

xmlfile='Glycolysis';
xlsname='Glycolysis_Constraints.xlsx';

%% create an EMU model
tic;
[model,memu]=modelinit(xmlfile,{'C','H'});
simulatefcn=str2func(xmlfile);
toc;

%% read constraints

free_net=zeros(size(model.kernel_net,2),1);
free_xch=zeros(size(model.kernel_xch,2),1);

[ineq,ineqStr]=xlsreadineq(xlsname,model,free_net,free_xch);
[eq,eqStr]=xlsreadeq(xlsname,model,free_net,free_xch);

%% randomize inital fluxes within constraints
n=100000; % the number of sets of initial fluxes for each input metabolite labeling
seed=0;
[free_net_set,free_xch_set]=fluxinit(model,ineq,eq,n,1,seed);
free_set=[free_net_set;free_xch_set]'; % OUTPUT for ML

%% simulate isotope labeling 
tracers = [1,1,0,0,0,0,0,0,0,0,1,0,0];                  % List out tracer configs used
for j1=1:size(tracers,1)                                % iterate through all tracers
    input.GLC_IN__C=isotopomervector(tracers(j1,1:6),1);
    input.GLC_IN__H=isotopomervector(tracers(j1,:),1);
    input.H_IN__H=isotopomervector([0],1);
    input.NADH_IN__H=isotopomervector([0],1);
    clear(char(simulatefcn))                            % reset the input metabolites, which is persistent in fcn

    for i=n:-1:1                                        % run simulation with input metabolite
        simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
    end
    simulated_sets(j1)=mergestruct(simulated);          % INPUT for ML (one tracer)
end
[~,simulated_set,lmid]=mergestruct2(simulated_sets);         % INPUT for ML (multiple tracers in struct)
singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)
