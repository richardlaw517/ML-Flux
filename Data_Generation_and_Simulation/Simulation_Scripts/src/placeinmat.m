function [mat,vec,keepmet,keepexp] = placeinmat(tracer,mea,exp,alltracers,mets,lmid,nanor)
%% placeinmat places user measurement (isotope tracer experiments and metabolite labeling) into correct places in data matrix
% input: tracer, a matrix listing the tracers used for isotope tracing
% experiments; mea, a struct showing measured metabolites and labeling
% fractions; exp, names of experiments; alltracers, the matrix listing all 
% of the tracers that are used for the training ML model; nmet, the number 
% of metabolites; lmid, the length of the MID vector of the largest 
% metabolite; and nanor, a placeholder for removed data (default NaN).
% output: a matrix and a vector whose elements show measured data in
% correct positions and NaN (or another placeholder) in the rest of the
% matrix and the vector as well as their positions.

if nargin<6 || isempty(nanor)
    nanor=NaN;
end

nexp=size(tracer,1);
nmet=length(mets);
keepmet=zeros(nexp,nmet);
mat=nanor*ones(size(alltracers,1)*lmid,nmet);

for i=nexp:-1:1
    keepexp(i)=find(ismember(alltracers,tracer(i,:),'rows'));
    mets_mea=fieldnames(mea);
    mets_temp=mets;
    expind=cell2mat(strfind(mets,'__'))+2;
    for j=1:nmet
        mets_temp{j}(expind(j):end)='';
        mets_temp{j}=[mets_temp{j} exp{i}];
    end
    for k=1:length(mets_mea)
        metpos=find(ismember(mets_temp,mets_mea{k}));
        if ~isempty(metpos)
            keepmet(i,metpos)=1;
            mat((keepexp(i)-1)*lmid+1:keepexp(i)*lmid,metpos)=mea.(mets_mea{k});            
        end
    end
    vec=reshape(mat,1,[]);
end