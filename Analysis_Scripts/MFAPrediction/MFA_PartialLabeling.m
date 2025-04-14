%% Model setup

model=load("CCM.mat").mod;
model.mets{'Suc'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));
model.mets{'Fum'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));
for i= 1:length(model.expts)
    model.expts(i).data_ms(2).id='GLC';
    model.expts(i).data_ms(22).id='AcCOA';
    model.expts(i).data_ms(25).id='SuccCOA';
end
TransportReactions={'Glucose_Upt','CO2_Upt','CO2_Exp','Gln_Ex','LDH','OAA_Upt','Glu_Upt','Gln_Upt','OGA_Glu','AcCoA_Ex','DHAP_Exp','PGA_Exp','PYR_Exp','LAC_Exp','G6P_Exp','R5P_Ex','OAA_Ex','CoA_Uptk','ALA_synth','ALA_Ex'};
model.rates.flx.fix=0;
model.rates{'Glucose_Upt'}.flx.fix = 1;
names_order=model.rates{'Glucose_Upt','CO2_Upt','CO2_Exp','Gln_Ex','LDH','OAA_Upt','Gln_Upt','Glu_Upt','OGA_Glu','AcCoA_Ex','DHAP_Exp','PGA_Exp','PYR_Exp','LAC_Exp','G6P_Exp','R5P_Ex','OAA_Ex','CoA_Uptk','ALA_synth','ALA_Ex'}.id;
dict=dictionary(TransportReactions,1:20);
reordered_ind=dict(names_order);

data_transport=readtable("CCM_TransportFluxes.dat");

%% Fix specific transport fluxes and set constraints 
model.rates{'Glucose_Upt'}.flx.fix = 1;
model.rates{'CO2_Upt'}.flx.lb = 490;
model.rates{'CO2_Upt'}.flx.ub = 500;
model.rates{'CO2_Exp'}.flx.lb = 490;
model.rates{'CO2_Exp'}.flx.ub = 500;
model.rates{'Gln_Upt'}.flx.ub = 500;

%% Parallelized flux estimation with missing data

fluxEstimations=cell([11708 4]);
tic
for batchNum = 1001:11708
    opts = detectImportOptions("label_test_CCM.dat"); % Detect import options
    opts.DataLines = [batchNum, batchNum]; % Specify row range to read
    labelData = readtable("label_test_CCM.dat", opts); % Read with the specified options

    tracerExperiments=model.expts;
    %For each column get num atoms for this
    for i= 1:length(model.expts(1).data_ms)
        %Going through each metabolite
        for j=1:length(tracerExperiments)
        %Going through each tracer
            a=length(model.mets{(char(tracerExperiments(j).data_ms(i).id))}.atoms);
            temp=idv(table2array(labelData(1:end,1+ 8*20*(i-1) + 8*(j-1):1+ 8*20*(i-1) + 8*(j-1)+a))');
            temp.std=ones(size(temp.val))*1e-3;
            tracerExperiments(j).data_ms(i).idvs=temp;
        end
    end
    
    %Run the flux estimation
    modelupd=model;
    modelupd.expts=tracerExperiments;
    modelupd.expts.data_ms.idvs.on=~modelupd.expts.data_ms.idvs.on;
    for j=1:length(tracerExperiments(1).data_ms(1).idvs)
        %for each j we just turn on one of the experiments at the jth index
        Tracers_n=ceil(rand()*1); %Generates a distribution of numbers from 1 to 3
        Tracer_ind=randperm(length(modelupd.expts),Tracers_n); %Selects a random permutation of Tracer_n numbers in the range 
        Metabolites_n=ceil((rand()*21))+9; %Generates a distribution of numbers from 10 to 30
        Metabolites_ind=randperm(length(tracerExperiments(1).data_ms),Metabolites_n); %Selects a random permutation of Metabolites_n numbers in the range  
        temparr=(modelupd.expts.data_ms.idvs.on)*0;
        for i=Tracer_ind %Loop through only the selected tracer experiments
            for k=Metabolites_ind %Loop through only selected metabolites for each experiment
                index=(i-1)*length(tracerExperiments(1).data_ms.idvs)+(k-1)*length(tracerExperiments(1).data_ms(1).idvs);
                temparr(index+j)=1;
            end
        end
        modelupd.expts.data_ms.idvs.on=logical(temparr);
        transport_flux=data_transport(batchNum,:);
        %Recreate an experiment while adding the flux information
        F=data('Glucose_Upt CO2_Upt CO2_Exp Gln_Ex LDH OAA_Upt Glu_Upt Gln_Upt OGA_Glu AcCoA_Ex DHAP_Exp PGA_Exp PYR_Exp LAC_Exp G6P_Exp R5P_Ex OAA_Ex CoA_Uptk ALA_synth ALA_Ex');
        F.val=table2array(transport_flux);
        for i=1:length(tracerExperiments)
            modelupd.expts(i).data_flx=F;
        end
        fit=estimate(modelupd);

    end
    fluxEstimations(batchNum,:) = {fit.par.val(1:86),fit.par.std(1:86),tracerExperiments.id(Tracer_ind),tracerExperiments(1).data_ms.id(Metabolites_ind)};
end
toc
save("MFA_Predictions.mat")

