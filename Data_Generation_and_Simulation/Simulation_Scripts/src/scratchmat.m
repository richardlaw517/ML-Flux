function [scmat,keepmet,keepexp,mask] = scratchmat(mat,nmet,lmid,seed,nexp,nanor,bias)
%% scratchmat removes isotope tracer experiments and unmeasured metabolites
% input: a matrix (or vector) containing metabolite isotope labeling data
% from the full isotope tracer set; nmet, the number of metabolites; lmid, 
% the length of the MID vector of the largest metabolite; random number 
% seed (default 0); nexpt, the number of isotope tracer experiments to keep
% (default 1); nanor, a place holder for removed data (default NaN); and 
% bias, which will increase or decrease the number of measured metabolites
% depending on its value being greater than or less than 1 (default 1).
% output: a matrix (or vector) whose elements have been randomly replaced
% by NaN (or nanor) in such a way that mimics real situations that have one
% (or few) isotope tracer(s) has been used and a fraction of metabolites 
% are measured; keepmet, metabolites that are kept; keepexp, experiments 
% that are kept; and mask, a matrix of 0s and 1s that represent kept pixels
% and removed pixels, respectively.

if nargin<4 || isempty(seed)
    seed=0;
end
if nargin<5 || isempty(nexp)
    nexp=1;
end
if nargin<6 || isempty(nanor)
    nanor=NaN;
end
if nargin<7 || isempty(bias)
    bias=1;
end
rng(seed);

niso=size(mat,2)/nmet/lmid;
nsim=size(mat,1);
keepexp=randi(niso,nsim,nexp);
keepmet=bias*round(rand(nsim,nmet));

for i=nsim:-1:1
    temp=reshape(mat(i,:),[],nmet);
    scvec=nanor*ones(size(temp));
    if sum(keepmet(i,:))==0
        keepmet(i,randi(nmet))=1;
    end
    for j=1:nexp
        while max(keepexp(i,j)==keepexp(i,1:j-1))
            keepexp(i,j)=randi(niso);
        end
        for k=1:nmet
            if keepmet(i,k)>0
                scvec((keepexp(i,j)-1)*lmid+1:keepexp(i,j)*lmid,k)=temp((keepexp(i,j)-1)*lmid+1:keepexp(i,j)*lmid,k);
            end
        end
    end
    scmat(i,:)=reshape(scvec,1,[]);
    mask(i,:)=zeros(size(scmat(i,:)));
    mask(i,scmat(i,:)==nanor)=1;
end