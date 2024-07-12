%% labeling simulation script
xmlfile='GlyPPP';
xlsname='GlyPPP_Constraints.xlsx';

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
total_fluxes = 1000000;
n = 1000;
seed = 0;
free_net_set = zeros(size(model.kernel_net,2),total_fluxes);
free_xch_set = zeros(size(model.kernel_xch,2),total_fluxes);
for iterPar = 1:total_fluxes/n
    [free_net_set(1:end,n*(iterPar-1)+1:n*iterPar),free_xch_set(1:end,n*(iterPar-1)+1:n*iterPar)]=par_fluxinit(model,ineq,eq,n,1,seed);
    disp(strcat(['Finished ',num2str(iterPar),' out of ',num2str(1000)]))
end
free_set=[free_net_set;free_xch_set]'; % OUTPUT for ML

%% simulate isotope labeling 
input.CO2_IN__U=isotopomervector([0],1);

tracers=[1,1,1,0,0,0,0,0,0,0,0,0,0;
         1,1,0,0,0,0,0,0,0,0,0,0,0;
         1,0,0,0,0,0,0,0,0,0,0,0,0;
         0,1,1,1,1,1,0,0,0,0,0,0,0;
         0,1,1,0,0,0,0,0,0,0,0,0,0;
         0,1,0,0,0,0,0,0,0,0,0,0,0;
         0,0,1,1,0,0,0,0,0,0,0,0,0;
         0,0,1,0,0,0,0,0,0,0,0,0,0;
         0,0,0,1,1,1,0,0,0,0,0,0,0;
         0,0,0,1,1,0,0,0,0,0,0,0,0;
         0,0,0,1,0,0,0,0,0,0,0,0,0;
         0,0,0,0,1,0,0,0,0,0,0,0,0;
         0,0,0,0,0,1,0,0,0,0,0,0,0;
         1,1,1,1,1,1,0,0,0,0,0,0,0;
         0,1,0,0,1,1,0,0,0,0,0,0,0;
         1,0,0,0,1,1,0,0,0,0,0,0,0;
         1,1,0,1,0,0,0,0,0,0,0,0,0;
         0,1,0,1,0,0,0,0,0,0,0,0,0;
         1,0,0,1,0,0,0,0,0,0,0,0,0;
         0,0,0,1,0,1,0,0,0,0,0,0,0;
         0,1,0,0,1,0,0,0,0,0,0,0,0;
         1,0,1,0,0,0,0,0,0,0,0,0,0;
         0,0,0,0,1,1,0,0,0,0,0,0,0;
         1,0,0,0,0,1,0,0,0,0,0,0,0;];%List out tracer configs used
                                         
for tracerConfig = 1:size(tracers,1)                       
    tic
    input.GLC_IN__U=isotopomervector(tracers(tracerConfig,:),1);
    clear(char(simulatefcn))    % reset the input metabolite, which is persistent in fcn

    parfor i=1:total_fluxes                                        % run simulation with input metabolite
        simulated(i)=simulatefcn(free_net_set(:,i),free_xch_set(:,i),input);
    end

    simulated_sets(tracerConfig)=mergestruct(simulated);          % INPUT for ML (one tracer)
    disp(strcat(['Finished ',num2str(tracerConfig),' fluxes out of ',num2str(size(tracers,1))]))
    toc
end

[~,simulated_set,lmid]=mergestruct2(simulated_sets); 
singlemat=struct2mat(simulated_set);                    % INPUT for ML (multiple tracers in a single matrix)