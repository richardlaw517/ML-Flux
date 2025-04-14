%% Model setup

model=load("CCM.mat").mod;

model.mets{'Suc'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));
model.mets{'Fum'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));
for i= 1:length(model.expts)
   model.expts(i).data_ms(2).id='GLC';
   model.expts(i).data_ms(22).id='AcCOA';
   model.expts(i).data_ms(25).id='SuccCOA';
end

results2=cell([11708 2]);

for batchNum = 1:1
    tic
    opts = detectImportOptions("label_test_CCM.dat"); % Detect import options
    opts.DataLines = [batchNum, batchNum]; % Specify row range to read
    data2 = readtable("label_test_CCM.dat", opts); % Read with the specified options

    Experiments2=model.expts;

    for i= 1:length(model.expts(1).data_ms)
        %Going through each metabolite
        for j=1:length(Experiments2)
        %Going through each tracer
            % k = 1+ 8*20*(i-1) + 8*(j-1)
            a=length(model.mets{(char(Experiments2(j).data_ms(i).id))}.atoms); %Test it once
            % kupd=k+a;
            temp=idv(table2array(data2(1:end,1+ 8*20*(i-1) + 8*(j-1):1+ 8*20*(i-1) + 8*(j-1)+a))');
            temp.std=ones(size(temp.val))*1e-3;
            Experiments2(j).data_ms(i).idvs=temp;
            % k=k+8;
        end
    end
        
    modelupd=model;
    modelupd.expts=Experiments2;
    modelupd.expts.data_ms.idvs.on=~modelupd.expts.data_ms.idvs.on;
    for j=1:length(Experiments2(1).data_ms(1).idvs)
        %for each j we just turn on one of the experiments at the jth index
        temparr=(modelupd.expts.data_ms.idvs.on)*0;
        for i=1:length(Experiments2)
            for k=1:length(Experiments2(1).data_ms)
                index=(i-1)*length(Experiments2(1).data_ms.idvs)+(k-1)*length(Experiments2(1).data_ms(1).idvs);
                temparr(index+j)=1;
            end
        end
        modelupd.expts.data_ms.idvs.on=logical(temparr)
        fit=estimate(modelupd);
    end
    results2(batchNum,:) = {fit.par.val(1:86),fit.par.std(1:86)};
    toc
end