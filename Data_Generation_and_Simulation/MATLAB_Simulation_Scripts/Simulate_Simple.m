%% labeling simulation script
xmlfile='Simple';
xlsname='Simple_Constraints.xlsx';

%% create an EMU model
tic;
[model,memu]=modelinit(xmlfile,{'U'});
simulatefcn=str2func(xmlfile);
toc;

%% read constraints
free_net=zeros(size(model.kernel_net,2),1);
free_xch=zeros(size(model.kernel_xch,2),1);

[ineq,ineqStr]=xlsreadineq(xlsname,model,free_net,free_xch);
[eq,eqStr]=xlsreadeq(xlsname,model,free_net,free_xch);

%% randomize inital fluxes within constraints
n=100; % the number of sets of initial fluxes for each input metabolite labeling
seed=0;
[free_net_set,free_xch_set]=fluxinit_logUniXch(model,ineq,eq,n,1,seed);
free_set=[free_net_set;free_xch_set]'; % OUTPUT for ML

%% simulate isotope labeling 
m1=3;                                                   % # of atoms in input metabolite 1; may need multiple inpute mets
for j1=2^m1-2:-1:1                                      % generate all nontrivial tracers
    alltracers(j1,:)=de2bi(j1,m1);
    input.molA__U=isotopomervector(alltracers(j1,:),1);
    clear(char(simulatefcn))                            % reset the input metabolite, which is persistent in fcn
    for i=n:-1:1                                        % run simulation with input metabolite
        simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
    end
    simulated_sets(j1)=mergestruct(simulated);          % INPUT for ML (one tracer)
end
[~,simulated_set,lmid]=mergestruct2(simulated_sets);    % INPUT for ML (multiple tracers in struct)
singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)