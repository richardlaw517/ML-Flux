function [m,mexp]=modelinit(file,exp,metfrag,optfilename)
%% Initiate model
% file is .xml file with metabolite and reaction definition
% metfrag is metabolites and fragments to simulate and fit

% try
%% Read XML file to create stoichiometric model
if isempty(strfind(file,'.'))
	file=[file '.xml'];
end
if nargin<4 || isempty(optfilename)
	if isempty(strfind(file,'.'))
		optfilename=file;
	else
		optfilename=file(1:strfind(file,'.')-1);
	end
end

if nargin<2
	exp=[];
elseif isnumeric(exp)
	exp=num2cell(exp);
end
exptype=zeros(size(exp));
atoms={'C' 'H' 'O' 'N'}; % C=1 H=2 O=4 N=8
for i=1:length(exp)
	exp{i}=num2str(exp{i});
	for a=1:length(atoms)
		if ~isempty(strfind(exp{i},atoms{a}))
			exptype(i)=exptype(i)+2^(a-1);
		end
	end
	if exptype(i)==0
		exptype(i)=1;
	end
end
[uexptype,~,exptypei]=unique(exptype);
% % % check if file exists already and stop if so
% if isempty(strfind(optfilename,'.'))
% 	if exist([optfilename '.m'],'file')
% 		errordlg(['File already exists: ' [optfilename '.m']],'File Error');
% % 		warningMessage = sprintf('Warning: file already exists:\n%s',[optfilename '.m']);
% % 		uiwait(msgbox(warningMessage));
% 		m=[];
% 		return
% 	end
% else
% 	if exist(optfilename,'file')
% 		errordlg(['File already exists: ' optfilename],'File Error');
% % 		warningMessage = sprintf('Warning: file already exists:\n%s',optfilename);
% % 		uiwait(msgbox(warningMessage));
% 		m=[];
% 		return
% 	end
% end

xml=xml2struct2(file);
xml=xml.fluxml;

% % Parse metabolites
metss=xml.reactionnetwork.metabolitepools.pool;
metss=cellfun(@(x) x.Attributes,metss,'UniformOutput',false);
for i=1:length(metss)
	if isfield(metss{i},'atoms')
		metss{i}.atoms=str2double(metss{i}.atoms);
		metss{i}.C=metss{i}.atoms;
	elseif isfield(metss{i},'C')
		metss{i}.C=str2double(metss{i}.C);
		metss{i}.atoms=metss{i}.C;
    else
        metss{i}.C=0;
		metss{i}.atoms=0;
	end
	for j=2:length(atoms)
		if isfield(metss{i},atoms{j})
			metss{i}.(atoms{j})=str2double(metss{i}.(atoms{j}));
			metss{i}.atoms=metss{i}.atoms+metss{i}.(atoms{j});
		else
			metss{i}.(atoms{j})=0;
		end	
	end
	if isfield(metss{i},'symm') % angle of symmetry: 360/symm = # of equivs
		metss{i}.symm=str2double(metss{i}.symm);
	else
		metss{i}.symm=360;
	end
end

% % Parse reactions
rxnss=xml.reactionnetwork.reaction;
S=zeros(length(metss),length(rxnss));
rxn=ones(length(rxnss),1);%bug fix JP4 10/5/21
minS=S; % used for reaction with a metabolite that is both a educt and a product
maxS=S;
equiv=zeros(length(rxnss),1);
leak=zeros(length(rxnss),1);

for i=1:length(rxnss)
	rxnss{i}.id=rxnss{i}.Attributes.id;
	if isfield(rxnss{i}.Attributes,'equiv')
		rxnss{i}.equiv=strsplit(rxnss{i}.Attributes.equiv);
		equiv(i)=length(rxnss{i}.equiv);
	end
	if isfield(rxnss{i}.Attributes,'leak')
		rxnss{i}.leak=strsplit(rxnss{i}.Attributes.leak);
		leak(i)=length(rxnss{i}.leak);
	end
	nr=0;
	np=0;
	if length(rxnss{i}.reduct)>1
		for j=1:length(rxnss{i}.reduct)
			rxnss{i}.reduct{j}=rxnss{i}.reduct{j}.Attributes;
			rxnss{i}.reduct{j}.index=find(cellfun(@(x) strcmp(x.id,rxnss{i}.reduct{j}.id),metss));
            if isempty(rxnss{i}.reduct{j}.index)
				m=[];
				mexp=[];
                error([rxnss{i}.reduct{j}.id ' reactant in reaction ' rxnss{i}.id ' is not defined.'])
				return
            end
			S(rxnss{i}.reduct{j}.index,i)=S(rxnss{i}.reduct{j}.index,i)-1;
			minS(rxnss{i}.reduct{j}.index,i)=minS(rxnss{i}.reduct{j}.index,i)-1;
			if length(rxnss{i}.reduct{j}.cfg)~=metss{rxnss{i}.reduct{j}.index}.atoms
				m=[];
				mexp=[];
                error([rxnss{i}.reduct{j}.id ' in reaction ' rxnss{i}.id ' has an error.'])
				return
			else
				nr=nr+length(rxnss{i}.reduct{j}.cfg);
			end
		end
	elseif length(rxnss{i}.reduct)==1
		rxnss{i}.reduct=rxnss{i}.reduct.Attributes;
		rxnss{i}.reduct.index=find(cellfun(@(x) strcmp(x.id,rxnss{i}.reduct.id),metss));
        if isempty(rxnss{i}.reduct.index)
            m=[];
            mexp=[];
            delete(h)
            error([rxnss{i}.reduct.id ' reactant in reaction ' rxnss{i}.id ' is not defined.'])
            return
        end
		S(rxnss{i}.reduct.index,i)=S(rxnss{i}.reduct.index,i)-1;
		minS(rxnss{i}.reduct.index,i)=minS(rxnss{i}.reduct.index,i)-1;
		if length(rxnss{i}.reduct.cfg)~=metss{rxnss{i}.reduct.index}.atoms
			m=[];
			mexp=[];
            error([rxnss{i}.reduct.id ' in reaction ' rxnss{i}.id ' has an error.'])
			return
		else
			nr=nr+length(rxnss{i}.reduct.cfg);
		end
	end
% 	minS(S<0)=S(S<0);
% 	maxS(S>0)=S(S>0);
	if isfield(rxnss{i},'rproduct')
		if length(rxnss{i}.rproduct)>1
			for j=1:length(rxnss{i}.rproduct)
				rxnss{i}.rproduct{j}=rxnss{i}.rproduct{j}.Attributes;
				rxnss{i}.rproduct{j}.index=find(cellfun(@(x) strcmp(x.id,rxnss{i}.rproduct{j}.id),metss));
                if isempty(rxnss{i}.rproduct{j}.index)
                    m=[];
                    mexp=[];
                    error([rxnss{i}.rproduct{j}.id ' product in reaction ' rxnss{i}.id ' is not defined.'])
                    return
                end
				S(rxnss{i}.rproduct{j}.index,i)=S(rxnss{i}.rproduct{j}.index,i)+1;
				maxS(rxnss{i}.rproduct{j}.index,i)=maxS(rxnss{i}.rproduct{j}.index,i)+1;
				if length(rxnss{i}.rproduct{j}.cfg)~=metss{rxnss{i}.rproduct{j}.index}.atoms
					m=[];
					mexp=[];
                    error([rxnss{i}.rproduct{j}.id ' in reaction ' rxnss{i}.id ' has an error.'])
					return
				else
					np=np+length(rxnss{i}.rproduct{j}.cfg);
				end
			end
		elseif length(rxnss{i}.rproduct)==1
			rxnss{i}.rproduct=rxnss{i}.rproduct.Attributes;
			rxnss{i}.rproduct.index=find(cellfun(@(x) strcmp(x.id,rxnss{i}.rproduct.id),metss));
            if isempty(rxnss{i}.rproduct.index)
                m=[];
                mexp=[];
                error([rxnss{i}.rproduct.id ' product in reaction ' rxnss{i}.id ' is not defined.'])
                return
            end
			S(rxnss{i}.rproduct.index,i)=S(rxnss{i}.rproduct.index,i)+1;
			maxS(rxnss{i}.rproduct.index,i)=maxS(rxnss{i}.rproduct.index,i)+1;
			if length(rxnss{i}.rproduct.cfg)~=metss{rxnss{i}.rproduct.index}.atoms
				m=[];
				mexp=[];
                error([rxnss{i}.rproduct.id ' in reaction ' rxnss{i}.id ' has an error.'])
				return
			else
				np=np+length(rxnss{i}.rproduct.cfg);
			end
		end
		if nr~=np
			m=[];
			mexp=[];
            error(['The number of atoms in reactants and products in reaction ' rxnss{i}.id ' are different.'])
			return
		end
	end
% 	minS(S<0)=S(S<0);
% 	maxS(S>0)=S(S>0);
	rxnss{i}=rmfield(rxnss{i},'Attributes');
end
minS(S>0)=S(S>0);
maxS(S<0)=S(S<0);

% % Delineate network properties
noninput=logical(max(S,[],2)+(sum(S,2)~=-1));
m.S_=S;
m.minS=minS;
m.maxS=maxS;
m.info=xml.info;
m.types_fluxes={'net','xch','DEPD.n','DEPD.x','QCON.n','QCON.x','FREE.n','FREE.x','CONS.n','CONS.x'};
m.metss=metss;
m.rxnss=rxnss;
m.equiv=logical(equiv);
m.leak=logical(leak);

S=S(noninput,:);
S=[S zeros(size(S,1),sum(leak))];
m.S=S;

mets=cellfun(@(x) x.id,m.metss,'UniformOutput',false);
m.mets_=mets;
m.mets=mets(noninput);
m.input=mets(~noninput);

rxns=cellfun(@(x) x.id,m.rxnss,'UniformOutput',false);
m.rxns_=rxns;
m.rxns=[rxns cellfun(@(x) [x.id '_leak'],rxnss(m.leak),'UniformOutput',false)];

fluxes=zeros(length(m.rxns),10);
for i=1:length(m.rxns)
	if min(S(:,i))<0 && max(S(:,i))>0
		fluxes(i,8)=1;
	elseif min(S(:,i))==0 && max(S(:,i))==0
		fluxes(i,8)=1;
	else
		fluxes(i,6)=1;
	end
end
kernel=eye(size(m.rxns,2));
m.kernel_xch=kernel(:,logical(fluxes(:,8)));

kernel=null(S,'r');
m.kernel_net=kernel(:,1:end-sum(leak));

r=rref(S);

for i=1:size(r,1)
	for j=1:size(r,2)
		if r(i,j)==1
			fluxes(j,3)=1; % dependent fluxes
			break
		end
	end
end
fluxes(end-sum(leak)+1:end,3)=1;
fluxes(:,7)=~fluxes(:,3);
m.fluxes=fluxes;


%% Search for EMU
if nargin<3 || isempty(metfrag)
	metfrag=m.mets;
elseif ~isstruct(metfrag)
%         metfrag=metfrag;
    if ~iscell(metfrag)
        % 	metfrag=[m.mets metfrag];
        metfrag={metfrag};
    end
end

% % Create EMU for each experiment
mexp=cell(size(uexptype));
if isempty(strfind(optfilename,'.'))
	fid_=fopen([optfilename '.m'],'w');
	fprintf(fid_,'function mol=%s(free_net,free_xch,inp)\n\n',optfilename);
else
	fid_=fopen(optfilename,'w');
	fprintf(fid_,'function mol=%s(free_net,free_xch,inp)\n\n',optfilename(1:strfind(optfilename,'.')-1));
end

% % wrapper function to simulate multiple labeling experiments
fprintf(fid_,'mol.EMU=1;\n');

for mexpi=1:length(uexptype)
mexp{mexpi}.id=exp(exptypei==mexpi);
mexp{mexpi}.exptype=cell2mat(atoms(logical(de2bi(uexptype(mexpi),length(atoms)))));

for i=length(metfrag):-1:1
	if ~isfield(metfrag{i},'cfg') % atom numeral indices
		metfrag_{i}.cfg=[];
	else
		metfrag_{i}.cfg=metfrag{i}.cfg;
    end
    if ~isfield(metfrag{i},'id')
		metfrag_{i}.id=metfrag{i};
        frag=strfind(metfrag_{i}.id,'_');
        if ~isempty(frag)
            frag=[frag length(metfrag_{i}.id)+1];
            metfrag_{i}.id=metfrag_{i}.id(1:frag(1)-1);
            if isempty(metfrag_{i}.cfg)
                for j=length(frag)-1:-1:1
                    metfrag_{i}.cfg(j)=str2double(metfrag{i}(frag(j)+1:frag(j+1)-1));
                end
            end
        end
	else
		metfrag_{i}.id=metfrag{i}.id;
	end
	if ~isfield(metfrag{i},'meti')
		metfrag_{i}.meti=find(strcmp(mets,metfrag_{i}.id));
	else
		metfrag_{i}.meti=metfrag{i}.meti;
	end
	if ~isfield(metfrag{i},'equiv')
		metfrag_{i}.equiv=0;
	else
		metfrag_{i}.equiv=metfrag{i}.equiv;
	end
	if ~isfield(metfrag{i},'network') % EMU network #
		metfrag_{i}.network=[];
	else
		metfrag_{i}.network=metfrag{i}.network;
	end
	if ~isfield(metfrag{i},'index') % index in EMU matrix
		metfrag_{i}.index=[];
	else
		metfrag_{i}.index=metfrag{i}.index;
	end
	% atoms: number of atoms
	if isempty(metfrag_{i}.cfg)
		metfrag_{i}.atoms=0;
		for j=find(de2bi(uexptype(mexpi),length(atoms)))
			metfrag_{i}.atoms=metfrag_{i}.atoms+metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j});
			if j==1
				metfrag_{i}.cfg=1:metfrag_{i}.atoms;
			else
				otheratoms=0;
				for k=1:j-1
					otheratoms=otheratoms+metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j-k});
				end
				metfrag_{i}.cfg=[metfrag_{i}.cfg otheratoms+(1:metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j}))];
			end
		end
    else
        relevantatoms=false(size(metfrag_{i}.cfg));
		for j=find(de2bi(uexptype(mexpi),length(atoms)))
			if j==1
				relevantatoms=relevantatoms | ismember(metfrag_{i}.cfg,1:metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j}));
% 				metfrag_{i}.cfg(ismember(metfrag_{i}.cfg,1:metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j})))=[];
			else
				otheratoms=0;
				for k=1:j-1
					otheratoms=otheratoms+metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j-k});
				end
				relevantatoms=relevantatoms | ismember(metfrag_{i}.cfg,otheratoms+(1:metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j})));
% 				metfrag_{i}.cfg(ismember(metfrag_{i}.cfg,otheratoms+(1:metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j}))))=[];
			end
		end
		metfrag_{i}.cfg(~relevantatoms)=[];
		metfrag_{i}.atoms=length(metfrag_{i}.cfg);
	end
	if ~isfield(metfrag{i},'input') % index in EMU matrix
		metfrag_{i}.input=0;
	else
		metfrag_{i}.input=metfrag{i}.input;
	end
end
emusize=cellfun(@(x) x.atoms,metfrag_);
metfrag_=metfrag_(emusize~=0);
mexp{mexpi}.metfrag=metfrag_;

nemu=0;
i=max(emusize);

maxemu=max(1700,round(2^i*sqrt(length(rxns))*2^length(mexp{mexpi}.exptype)));
emu=cell((maxemu),1);
to=cell((maxemu),1);
from=cell((maxemu),1);
forward=cell((maxemu),1);
reverse=cell((maxemu),1);

add2emu=metfrag_;
while ~isempty(add2emu) || i>0
	last=nemu+length(add2emu);
	emu(nemu+1:last)=add2emu;
	emusize=cellfun(@(x) max(x.atoms),emu(1:last));
	cc=find(emusize==i);
	nemu=last;
	for j=cc'
		if isempty(to{j}) && isempty(from{j})
			id=emu{j}.id;
			cfg=emu{j}.cfg;
			symm=metss{cellfun(@(x) strcmp(x.id,id),metss)}.symm;
			nsymm=360/symm;
			C=metss{cellfun(@(x) strcmp(x.id,id),metss)}.C;
			H=metss{cellfun(@(x) strcmp(x.id,id),metss)}.H;
			O=metss{cellfun(@(x) strcmp(x.id,id),metss)}.O;
			N=metss{cellfun(@(x) strcmp(x.id,id),metss)}.N;
			tempR=maxS(strcmp(mets,id),:);
			forward{j}=find(tempR>0);
			tempR=minS(strcmp(mets,id),:);
			reverse{j}=find(tempR<0);
			addemu_1=[];
			addemu_2=[];
			addemu_3=[];
			lasto1=0;
			lasto2=0;
			lasto3=0;
			for k=forward{j} %reactions that produce the current met
				tempC=minS(:,k);
				l=find(tempC<0); %precursors to the current met
				from{j}=[from{j}; l];
				if length(rxnss{k}.rproduct)>1
					kk=find(cellfun(@(x) strcmp(x.id,id),rxnss{k}.rproduct));
				else
					kk=1;
				end
				for ll=kk
					if length(rxnss{k}.rproduct)>1
						curr=rxnss{k}.rproduct{ll}.cfg(cfg); %current met atom letters
						currfull=rxnss{k}.rproduct{ll}.cfg;
					else
						curr=rxnss{k}.rproduct.cfg(cfg); %current met atom letters
						currfull=rxnss{k}.rproduct.cfg;
					end
					if equiv(k)>0 || nsymm>1
						cfg_=cfg;
						for n=equiv(k):-1:1
							ecurrcfgl=rxnss{k}.equiv{n};
							ecurratm=find(ismember(currfull,ecurrcfgl));
							ecurrbool=ismember(ecurrcfgl,currfull);
							if ecurrbool
								cfg_(ismember(cfg_,ecurratm))=[]; % remove identical atm from cfg
								ecurrn=sum(ismember(curr,ecurrcfgl));
								ecurratm_{n}=combnk(ecurratm,ecurrn);
							else
								ecurratm_{n}=[];
							end
						end
						ecurratm=1;
						for n=1:equiv(k) % permutation for more than 1 equiv
							if ~isempty(ecurratm_{n})
								ecurratm=combvec(ecurratm,ecurratm_{n}');
							end
						end
						ecurratm(1,:)=[];
						ecurratm=ecurratm';
						p=size(ecurratm,1);
						for o=p:-1:1
							if ~isempty(ecurratm) && nsymm==1
								newcfg=sort([cfg_ ecurratm(o,:)]);
								equivj=find(cellfun(@(x) strcmp(x.id,id),emu(1:last)) & cellfun(@(x) isequal(x.cfg,newcfg),emu(1:last)));
								if isempty(equivj)  %check if equiv emus are already in emu. if not, add them
									addemu_1{lasto1+o}.id=id;
									addemu_1{lasto1+o}.cfg=newcfg;
									addemu_1{lasto1+o}.meti=emu{j}.meti;
									addemu_1{lasto1+o}.equiv=j;
									addemu_1{lasto1+o}.network=0;
									addemu_1{lasto1+o}.index=0;
									addemu_1{lasto1+o}.atoms=length(newcfg);
									addemu_1{lasto1+o}.input=emu{j}.input;
								elseif sum(ismember(currfull,newcfg))<length(currfull)
									emu{equivj}.equiv=j;
								end
							elseif ~isempty(ecurratm) && nsymm==2
								newcfg=sort([cfg_ ecurratm(o,:)]);
								newcfg(newcfg<=C)=C+1-newcfg(newcfg<=C);
								newcfg(newcfg>C & newcfg<=C+H)=C+C+H+1-newcfg(newcfg>C & newcfg<=C+H);
								newcfg(newcfg>C+H & newcfg<=C+H+O)=C+C+H+H+O+1-newcfg(newcfg>C+H & newcfg<=C+H+O);
								newcfg(newcfg>C+H+O & newcfg<=C+H+O+N)=C+C+H+H+O+O+N+1-newcfg(newcfg>C+H+O & newcfg<=C+H+O+N);
								newcfg=sort(newcfg);
								equivj=find(cellfun(@(x) strcmp(x.id,id),emu(1:last)) & cellfun(@(x) isequal(x.cfg,newcfg),emu(1:last)));
								if isempty(equivj)  %check if equiv emus are already in emu. if not, add them
									addemu_1{lasto1+p+o}.id=id;
									addemu_1{lasto1+p+o}.cfg=newcfg;
									addemu_1{lasto1+p+o}.meti=emu{j}.meti;
									addemu_1{lasto1+p+o}.equiv=j;
									addemu_1{lasto1+p+o}.network=0;
									addemu_1{lasto1+p+o}.index=0;
									addemu_1{lasto1+p+o}.atoms=length(newcfg);
									addemu_1{lasto1+p+o}.input=emu{j}.input;
								elseif sum(ismember(currfull,newcfg))<length(currfull)
									emu{equivj}.equiv=j;
								end
							elseif isempty(ecurratm) && nsymm==2
% 								newcfg=sort([cfg_ ecurratm(o,:)]);
								newcfg=cfg_;
								newcfg(newcfg<=C)=C+1-newcfg(newcfg<=C);
								newcfg(newcfg>C & newcfg<=C+H)=C+C+H+1-newcfg(newcfg>C & newcfg<=C+H);
								newcfg(newcfg>C+H & newcfg<=C+H+O)=C+C+H+H+O+1-newcfg(newcfg>C+H & newcfg<=C+H+O);
								newcfg(newcfg>C+H+O & newcfg<=C+H+O+N)=C+C+H+H+O+O+N+1-newcfg(newcfg>C+H+O & newcfg<=C+H+O+N);
								newcfg=sort(newcfg);
								equivj=find(cellfun(@(x) strcmp(x.id,id),emu(1:last)) & cellfun(@(x) isequal(x.cfg,newcfg),emu(1:last)));
								if isempty(equivj)  %check if equiv emus are already in emu. if not, add them
									addemu_1{lasto1+p+o}.id=id;
									addemu_1{lasto1+p+o}.cfg=newcfg;
									addemu_1{lasto1+p+o}.meti=emu{j}.meti;
									addemu_1{lasto1+p+o}.equiv=j;
									addemu_1{lasto1+p+o}.network=0;
									addemu_1{lasto1+p+o}.index=0;
									addemu_1{lasto1+p+o}.atoms=length(newcfg);
									addemu_1{lasto1+p+o}.input=emu{j}.input;
								end
							end
							if nsymm>2
								error('Symmetry angles less than 180 are not defined yet.')
							end
						end
						lasto1=lasto1+p*nsymm;
					end

					for o=-sum(tempC(l)):-1:1 %bug fix JP 8/6/21
						if length(rxnss{k}.reduct)>1
							newid=rxnss{k}.reduct{o}.id;
							prevfull=rxnss{k}.reduct{o}.cfg;
						else
							newid=rxnss{k}.reduct.id;
							prevfull=rxnss{k}.reduct.cfg;
						end
						newcfg=find(ismember(prevfull,curr));
						if sum(cellfun(@(x) strcmp(x.id,newid),emu(1:last)) & cellfun(@(x) isequal(x.cfg,newcfg),emu(1:last)))==0 && ~isempty(newcfg) %check if equiv emus are already in emu. if not, add them
							addemu_2{lasto2+o}.id=newid;
							addemu_2{lasto2+o}.cfg=newcfg;
							addemu_2{lasto2+o}.meti=find(strcmp(mets,newid));
							addemu_2{lasto2+o}.equiv=0;
							addemu_2{lasto2+o}.network=0;
							addemu_2{lasto2+o}.index=0;
							addemu_2{lasto2+o}.atoms=length(newcfg);
							addemu_2{lasto2+o}.input=sum(strcmp(m.input,newid));
						end
						if leak(k)>0
							lo=length(l);
							newcfg=find(ismember(prevfull,strrep(curr,rxnss{k}.leak{1},'')));
							if sum(cellfun(@(x) strcmp(x.id,newid),emu(1:last)) & cellfun(@(x) isequal(x.cfg,newcfg),emu(1:last)))==0 && ~isempty(newcfg) %check if equiv emus are already in emu. if not, add them
								addemu_2{lasto2+lo+o}.id=newid;
								addemu_2{lasto2+lo+o}.cfg=newcfg;
								addemu_2{lasto2+lo+o}.meti=find(strcmp(mets,newid));
								addemu_2{lasto2+lo+o}.equiv=0;
								addemu_2{lasto2+lo+o}.network=0;
								addemu_2{lasto2+lo+o}.index=0;
								addemu_2{lasto2+lo+o}.atoms=length(newcfg);
								addemu_2{lasto2+lo+o}.input=sum(strcmp(m.input,newid));
							end
						end
					end
					lasto2=lasto2+length(l)*(leak(k)+1);
				end
			end
			
			for k=reverse{j} %reactions that consume the current met
				tempC=maxS(:,k);
				l=find(tempC>0); %derivatives of the current met
				to{j}=[to{j}; l];
				if length(rxnss{k}.reduct)>1
					kk=find(cellfun(@(x) strcmp(x.id,id),rxnss{k}.reduct));
				else
					kk=1;
				end
				for ll=kk
					if length(rxnss{k}.reduct)>1
						curr=rxnss{k}.reduct{ll}.cfg(cfg); %current met atom letters
					else
						curr=rxnss{k}.reduct.cfg(cfg); %current met atom letters
					end
					for o=length(l):-1:1
						if length(rxnss{k}.rproduct)>1
							newid=rxnss{k}.rproduct{o}.id;
							postfull=rxnss{k}.rproduct{o}.cfg;
						else
							newid=rxnss{k}.rproduct.id;
							postfull=rxnss{k}.rproduct.cfg;
						end
						newcfg=find(ismember(postfull,curr));
						if sum(cellfun(@(x) strcmp(x.id,newid),emu(1:last)) & cellfun(@(x) isequal(x.cfg,newcfg),emu(1:last)))==0 && ~isempty(newcfg) %check if equiv emus are already in emu. if not, add them
							addemu_3{lasto3+o}.id=newid;
							addemu_3{lasto3+o}.cfg=newcfg;
							addemu_3{lasto3+o}.meti=find(strcmp(mets,newid));
							addemu_3{lasto3+o}.equiv=0;
							addemu_3{lasto3+o}.network=0;
							addemu_3{lasto3+o}.index=0;
							addemu_3{lasto3+o}.atoms=length(newcfg);
							addemu_3{lasto3+o}.input=sum(strcmp(m.input,newid));
						end
						if leak(k)>0
							lo=length(l);
							newcfg=find(ismember(postfull,strrep(curr,rxnss{k}.leak{1},'')));
							if sum(cellfun(@(x) strcmp(x.id,newid),emu(1:last)) & cellfun(@(x) isequal(x.cfg,newcfg),emu(1:last)))==0 && ~isempty(newcfg) %check if equiv emus are already in emu. if not, add them
								addemu_3{lasto3+lo+o}.id=newid;
								addemu_3{lasto3+lo+o}.cfg=newcfg;
								addemu_3{lasto3+lo+o}.meti=find(strcmp(mets,newid));
								addemu_3{lasto3+lo+o}.equiv=0;
								addemu_3{lasto3+lo+o}.network=0;
								addemu_3{lasto3+lo+o}.index=0;
								addemu_3{lasto3+lo+o}.atoms=length(newcfg);
								addemu_3{lasto3+lo+o}.input=sum(strcmp(m.input,newid));
							end
						end
					end
					lasto3=lasto3+length(l)*(leak(k)+1);
				end
			end
			add2emu=[addemu_1 addemu_2 addemu_3];
			if ~isempty(add2emu)
				add2emu=add2emu(~cellfun(@isempty,add2emu));
				[~,ia,ic]=unique(cellfun(@(x) x.id,add2emu,'UniformOutput',false),'stable');
				if length(ia)<length(ic)
					remove=false(size(add2emu));
					for iia=1:length(ia)
						if sum(ic==iia)>1
							for nia1=find(ic==iia)'
								for nia2=find(ic==iia)'
									if isequal(add2emu{nia1}.cfg,add2emu{nia2}.cfg) && nia1<nia2
										remove(nia2)=true;
									end
								end
							end
						end
					end
					add2emu(remove)=[];
				end
				break
			end
		else
			add2emu=[];
		end
	end
% 	emusize=cellfun(@(x) max(x.atoms),emu);
	if isempty(add2emu) %sum(emusize==i)==length(jj)
		i=i-1;
	end
end

emu=emu(1:last);
to=to(1:last);
from=from(1:last);
forward=forward(1:last);
reverse=reverse(1:last);

[lvl,I]=sort(cell2mat(cellfun(@(x) x.atoms, emu,'UniformOutput',false)));
emu=emu(I);
to=to(I);
from=from(I);
forward=forward(I);
reverse=reverse(I);

% m.equivemu=equivemu; 
mexp{mexpi}.to=to; % derivatives of the current emu
mexp{mexpi}.from=from; % precursors to the current emu
mexp{mexpi}.forward=forward; % reactions that produce the current emu
mexp{mexpi}.reverse=reverse; % reactions that consume the current emu

fault=cell2mat(cellfun(@(x) isempty(x), mexp{mexpi}.to,'UniformOutput',false)) & cell2mat(cellfun(@(x) isempty(x), mexp{mexpi}.from,'UniformOutput',false));
if sum(fault)>0
	error('Check metabolite definition or EMU may have been added by mistake')
	find(fault)
end


%% Create EMU network
l=length(lvl);
[lvl,I]=unique(lvl);
I=[I; l+1];
cl=40; %bug fix JP 9/7/21 %bug fix JP4 10/5/21
hemu=cellfun(@(x) isequal(x.id,'H'),emu);
% nnet=cl*length(lvl);

A=[];
B=[];
y=[];
yy=[];
for i=length(lvl):-1:1;
	llvl=I(i+1)-I(i);
	for j=cl:-1:1
		A{cl*(i-1)+j}=cell(llvl);
		B{cl*(i-1)+j}=cell(llvl,llvl*cl);
		y{cl*(i-1)+j}=cell(llvl*cl,1);
		yy{cl*(i-1)+j}=[];
	end
end
yyy=1;
% S=m.S_;

for k=1:length(lvl)
	kk=(k-1)*cl+1;
	rrr=0;
	for r=I(k):(I(k+1)-1) %emu index
        llvl=I(k+1)-I(k);%bug fix JP4 10/5/21
		if emu{r}.input==0
% 			emu{r}.network=(k-1)*cl+1;
			rr=r-I(k)+1;
			for f=forward{r} %index of rxn that makes current emu
				stoi=abs(maxS(emu{r}.meti,f));
                stoichiometry=stoi; % stoichiometry is the absolute value of stoichiometry %bug fix JP2 10/1/21
				rrr=rrr+1;
				Aleak=[];
				yAleak=[];
				if stoi==1
					A{kk}{rr,rr}=[A{kk}{rr,rr} '-b(' num2str(f) ')'];
				elseif stoi>1
					A{kk}{rr,rr}=[A{kk}{rr,rr} '-' num2str(stoi) '*b(' num2str(f) ')'];
				end
				prevf=zeros(I(k),4); % to remove same metabolite EMUs going into one row in y
				lll=1;
				for c=find(cell2mat(cellfun(@(x) sum(x.meti==from{r}),emu,'UniformOutput',false)))' %finding the indices of emu's that generate current emu 
					if emu{c}.input==1 || c>=I(k+1)
						continue
					end
					cc=c-I(k)+1;
% 					stoi=abs(S(emu{c}.meti,f));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					ro=r;
					if emu{r}.equiv==0
						stoi=1; % number of equivalent EMUs
						requiv=r;
					else
						requiv=find(cell2mat(cellfun(@(x) x.equiv==emu{r}.equiv,emu,'UniformOutput',false)))';
						stoi=length(requiv);
					end
					for r=requiv
						if length(rxnss{f}.rproduct)>1
							lp=find(cellfun(@(x) strcmp(x.id,emu{r}.id),rxnss{f}.rproduct)); %lp is the index(indices) of rproduct in the f(th) reaction that matches the current emu
						else
							lp=find(strcmp(rxnss{f}.rproduct.id,emu{r}.id));
						end
						if length(rxnss{f}.reduct)>1
							le=find(cellfun(@(x) strcmp(x.id,emu{c}.id),rxnss{f}.reduct));
						else
							le=find(strcmp(rxnss{f}.reduct.id,emu{c}.id));
	%						le=strfind(rxnss{f}.reduct.id,emu{c}.id);
						end
						for llp=lp
							for lle=le
                                stoi2=0;%bug fix JP4 10/5/21
								if length(rxnss{f}.reduct)>1 && length(rxnss{f}.rproduct)>1
									if ismember(rxnss{f}.reduct{lle}.cfg(emu{c}.cfg),rxnss{f}.rproduct{llp}.cfg(emu{r}.cfg))
										if cc<1
											prevf(lll,1)=c;
											prevf(lll,2)=emu{c}.meti;
											prevf(lll,3)=sum(ismember(rxnss{f}.reduct{lle}.cfg(emu{c}.cfg),rxnss{f}.rproduct{llp}.cfg(emu{r}.cfg)));%may need a fix in the future JP 8/6/21 % # of overlapping atoms
											prevf(lll,4)=r;
											lll=lll+1;
                                        else
                                            if stoichiometry>1 && length(rxnss{f}.rproduct{llp}.cfg)>=length(rxnss{f}.reduct{lle}.cfg(emu{c}.cfg)) && sum(cellfun(@(x) strcmp(x.id,rxnss{f}.rproduct{llp}.id),rxnss{f}.rproduct))==stoichiometry && sum(cellfun(@(x) strcmp(x.id,rxnss{f}.reduct{lle}.id),rxnss{f}.reduct))==stoichiometry%bug fix JP4 10/5/21
                                                stoi2=stoichiometry;%bug fix JP4 10/5/21
                                                rxn(f)=stoi2;%bug fix JP5 10/5/21
%                                             end%bug fix JP4 10/5/21 %bug fix JP6 10/11/21
                                            elseif stoichiometry>1 && length(rxnss{f}.reduct{lle}.cfg)==length(rxnss{f}.rproduct{llp}.cfg) && sum(cellfun(@(x) strcmp(x.id,rxnss{f}.rproduct{llp}.id),rxnss{f}.rproduct))==stoichiometry%bug fix JP6 10/11/21
                                                stoi2=stoichiometry;%bug fix JP6 10/11/21
                                            end%bug fix JP6 10/11/21
                                            if stoi==1
                                                if stoichiometry==1 || stoi2>0 || (stoichiometry>1 && length(rxnss{f}.rproduct{llp}.cfg)>length(rxnss{f}.reduct{lle}.cfg(emu{c}.cfg))) || (stoichiometry>1 && length(rxnss{f}.rproduct{llp}.cfg(emu{r}.cfg))<=length(rxnss{f}.reduct{lle}.cfg)/stoichiometry) %bug fix JP4 10/5/21 %bug fix JP8 10/19/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')'];%bug fix JP2 10/1/21
                                                else%bug fix JP2 10/1/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+' num2str(stoichiometry) '*f(' num2str(f) ')'];%bug fix JP2 10/1/21
                                                end%bug fix JP2 10/1/21
											elseif stoi>1
                                                if stoichiometry==1 || stoi2>0 || (stoichiometry>1 && length(rxnss{f}.rproduct{llp}.cfg)>length(rxnss{f}.reduct{lle}.cfg(emu{c}.cfg))) || (stoichiometry>1 && length(rxnss{f}.rproduct{llp}.cfg(emu{r}.cfg))<=length(rxnss{f}.reduct{lle}.cfg)/stoichiometry) %bug fix JP4 10/5/21 %bug fix JP8 10/19/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')/' num2str(stoi)];%bug fix JP2 10/1/21
                                                else%bug fix JP2 10/1/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')/' num2str(stoi/stoichiometry)];%bug fix JP2 10/1/21
                                                end%bug fix JP2 10/1/21
                                            end
										end
									end
								elseif length(rxnss{f}.reduct)>1 && length(rxnss{f}.rproduct)==1
									if ismember(rxnss{f}.reduct{lle}.cfg(emu{c}.cfg),rxnss{f}.rproduct.cfg(emu{r}.cfg))
										if cc<1
											prevf(lll,1)=c;
											prevf(lll,2)=emu{c}.meti;
											prevf(lll,3)=sum(ismember(rxnss{f}.reduct{lle}.cfg(emu{c}.cfg),rxnss{f}.rproduct.cfg(emu{r}.cfg)));
											prevf(lll,4)=r;
											lll=lll+1;
										else
											if stoi==1
												A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')'];
											elseif stoi>1
												A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')/' num2str(stoi)];
											end
										end
									end
								elseif length(rxnss{f}.reduct)==1 && length(rxnss{f}.rproduct)>1
									if [ismember(rxnss{f}.reduct.cfg(emu{c}.cfg),rxnss{f}.rproduct{llp}.cfg(emu{r}.cfg)) cc>0]
										if stoi==1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')'];
										elseif stoi>1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')/' num2str(stoi)];
										end
									end
								else
									if [ismember(rxnss{f}.reduct.cfg(emu{c}.cfg),rxnss{f}.rproduct.cfg(emu{r}.cfg)) cc>0]
										if stoi==1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')'];
										elseif stoi>1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+f(' num2str(f) ')/' num2str(stoi)];
										end
									end
								end
							end
						end
						if cc>0 && ~isempty(A{kk}{rr,cc})
							if A{kk}{rr,cc}(1)=='+'
								A{kk}{rr,cc}(1)=[];
							end
							if leak(f)>0 % only isomerase reaction for now
								leakcfg=find(ismember(rxnss{f}.reduct.cfg,strrep(rxnss{f}.reduct.cfg(emu{c}.cfg),rxnss{f}.leak{1},'')));
								if length(leakcfg)<lvl(k) && ~isempty(leakcfg)
% 									Aleak=A{kk}{rr,cc};
									Aleak=['-f(' num2str(f) ')'];
									A{kk}{rr,cc}=[A{kk}{rr,cc} '*(1-f(' num2str(length(leak)+sum(leak(1:f))) '))'];
									leakemu=cellfun(@(x) x.meti==emu{c}.meti & isequal(x.cfg,leakcfg),emu);
									yAleak=[',x' num2str(emu{leakemu}.network) '(' num2str(emu{leakemu}.index) '),x' num2str(emu{hemu}.network) '(' num2str(emu{hemu}.index) ')'];
								elseif isempty(leakcfg)
									hh=find(hemu)-I(k)+1;
									A{kk}{rr,cc}=[A{kk}{rr,cc} '*(1-f(' num2str(length(leak)+sum(leak(1:f))) '))'];
									A{kk}{rr,hh}=[A{kk}{rr,hh} '+f(' num2str(f) ')*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
									A{kk}{hh,cc}=[A{kk}{hh,cc} '+f(' num2str(f) ')*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
									A{kk}{hh,hh}=[A{kk}{hh,hh} '-f(' num2str(f) ')*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
% 									A{kk}{rr,cc}=[A{kk}{rr,cc} '*(1-f(' num2str(length(leak)+sum(leak(1:f))) '))'];
% 									yAleak=[',x' num2str(emu{hemu}.network) '(' num2str(emu{hemu}.index) ')'];
									if A{kk}{rr,hh}(1)=='+'
										A{kk}{rr,hh}(1)=[];
									end
									if A{kk}{hh,cc}(1)=='+'
										A{kk}{hh,cc}(1)=[];
									end
								end
							end
						end
					end
					r=ro;
				end
				if emu{r}.equiv==0
					stoi=1; % number of equivalent EMUs
					requiv=r;
				else
					requiv=find(cell2mat(cellfun(@(x) x.equiv==emu{r}.equiv,emu,'UniformOutput',false)))';
					stoi=length(requiv);
				end
				lvlf=zeros(stoi,1);
% 				Bleak=cell(stoi,1);
% 				yleak=cell(stoi,1);
				if sum(sum(prevf)) && isempty(cell2mat(regexp(A{kk}(rr,:),['f\(' num2str(f) '\)(?=$|[^/])'])))
					[uprevf,~,ic]=unique(prevf(:,[2 4]),'stable','rows');
                    [~,ia]=unique(prevf,'rows');                            %bug fix JP 8/6/21
                    iia=find(~ismember(1:max(ia),ia));                      %bug fix JP 8/6/21
                    rxnfstep=1;%bug fix JP5 10/5/21
                    if rxn(f)~=1 && ~isempty(cell2mat(regexp(A{kk}(rr,:),[num2str(rxn(f)) '\*b\(' num2str(f) '\)(?=$|[^/])'])))%bug fix JP5 10/5/21
                        rxnfstep=rxn(f);%bug fix JP5 10/5/21
                    end%bug fix JP5 10/5/21
                    if rxnfstep==1%bug fix JP5 10/5/21
                    for iiia=iia%bug fix JP 8/6/21
                        if prevf(iiia,3)>=lvl(k)-prevf(ia(end),3)%bug fix JP 8/12/21
                            rc=find(requiv==prevf(iiia,4));%bug fix JP 8/6/21
                            if prevf(iiia,3)<max(prevf(:,3)) %bug fix JP 8/6/21
                                y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(iiia)}.network) '(' num2str(emu{prevf(iiia)}.index) ')'];%bug fix JP 8/6/21
                                lvlf(rc)=lvlf(rc)+prevf(iiia,3);%bug fix JP 8/6/21
                            end%bug fix JP 8/6/21
                        end%bug fix JP 8/12/21
                    end%bug fix JP 8/6/21
                    end%bug fix JP4 10/5/21
					for iic=1:max(ic)
						if uprevf(iic,:)~=0
							mprevf=max(prevf(ic==iic,3));
							ci=find(ismember(prevf(:,2:4),[uprevf(iic,1) mprevf uprevf(iic,2)],'rows'));
							rc=find(requiv==uprevf(iic,2));
							for cii=1:rxnfstep:length(ci)%bug fix JP5 10/5/21
								if cii>1 && minS(emu{prevf(ci(cii))}.meti,f)<-1
									y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];
									lvlf(rc)=lvlf(rc)+mprevf;
								elseif cii==1
									y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];
									lvlf(rc)=lvlf(rc)+mprevf;
								end
							end
							if stoi==1
								B{kk}{rr,rrr+rc-1}=['-f(' num2str(f) ')'];
                                if rxnfstep>1%bug fix JP5 10/5/21
                                    B{kk}{rr,rrr+rc-1}=['-f(' num2str(f) ')*' num2str(rxn(f))];%bug fix JP5 10/5/21
                                end%bug fix JP4 10/5/21
							elseif stoi>1
								B{kk}{rr,rrr+rc-1}=['-f(' num2str(f) ')/' num2str(stoi)];
                                if rxnfstep>1%bug fix JP5 10/5/21
                                    B{kk}{rr,rrr+rc-1}=['-f(' num2str(f) ')/' num2str(stoi) '*' num2str(rxn(f))];%bug fix JP5 10/5/21
                                end%bug fix JP4 10/5/21
							end
% 							if leak(f)>0
% 								Bleak{rc}=B{kk}{rr,rrr+rc-1};
% 								B{kk}{rr,rrr+rc-1}=[B{kk}{rr,rrr+rc-1} '*(1-f(' num2str(length(leak)+sum(leak(1:f))) '))'];
% 							end
						end
					end
					for cii=1:stoi
						if lvlf(cii)~=lvl(k)
% 							Bleak{cii}=[];
% 							yleak{cii}=[];
							B{kk}{rr,rrr+cii-1}=[];
							y{kk}{rrr+cii-1}=[];
						end
					end
					rrr=rrr+stoi-1;
                end
                for cc=1:llvl%bug fix JP3 10/4/21
                    if sum(sum(prevf)) && ~isempty(cell2mat(regexp(A{kk}(rr,cc),['f\(' num2str(f) '\)(?=$|[^/])']))) && ~isempty(cell2mat(regexp(A{kk}(rr,cc),['b\(' num2str(f) '\)(?=$|[^/])'])))%bug fix JP3 10/4/21
                        [uprevf,~,ic]=unique(prevf(:,[2 4]),'stable','rows');%bug fix JP3 10/4/21
                        [~,ia]=unique(prevf,'rows');                            %bug fix JP3 10/4/21
                        iia=find(~ismember(1:max(ia),ia));                      %bug fix JP3 10/4/21
                        for iiia=iia%bug fix JP3 10/4/21
                            if prevf(iiia,3)>=lvl(k)-prevf(ia(end),3)%bug fix JP3 10/4/21
                                rc=find(requiv==prevf(iiia,4));%bug fix JP3 10/4/21
                                if prevf(iiia,3)<max(prevf(:,3)) %bug fix JP3 10/4/21
                                    y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(iiia)}.network) '(' num2str(emu{prevf(iiia)}.index) ')'];%bug fix JP3 10/4/21
                                    lvlf(rc)=lvlf(rc)+prevf(iiia,3);%bug fix JP3 10/4/21
                                end%bug fix JP3 10/4/21
                            end%bug fix JP3 10/4/21
                        end%bug fix JP3 10/4/21
                        for iic=1:max(ic)%bug fix JP3 10/4/21
                            if uprevf(iic,:)~=0%bug fix JP3 10/4/21
                                mprevf=max(prevf(ic==iic,3));%bug fix JP3 10/4/21
                                ci=find(ismember(prevf(:,2:4),[uprevf(iic,1) mprevf uprevf(iic,2)],'rows'));%bug fix JP3 10/4/21
                                rc=find(requiv==uprevf(iic,2));%bug fix JP3 10/4/21
                                for cii=1:length(ci)%bug fix JP3 10/4/21
                                    if cii>1 && minS(emu{prevf(ci(cii))}.meti,f)<-1%bug fix JP3 10/4/21
                                        y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];%bug fix JP3 10/4/21
                                        lvlf(rc)=lvlf(rc)+mprevf;%bug fix JP3 10/4/21
                                    elseif cii==1%bug fix JP3 10/4/21
                                        y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];%bug fix JP3 10/4/21
                                        lvlf(rc)=lvlf(rc)+mprevf;%bug fix JP3 10/4/21
                                    end%bug fix JP3 10/4/21
                                end%bug fix JP3 10/4/21
                                if stoi==1%bug fix JP3 10/4/21
                                    B{kk}{rr,rrr+rc-1}=['-f(' num2str(f) ')'];%bug fix JP3 10/4/21
                                elseif stoi>1%bug fix JP3 10/4/21
                                    B{kk}{rr,rrr+rc-1}=['-f(' num2str(f) ')/' num2str(stoi)];%bug fix JP3 10/4/21
                                end%bug fix JP3 10/4/21
                            end%bug fix JP3 10/4/21
                        end%bug fix JP3 10/4/21
                        for cii=1:stoi%bug fix JP3 10/4/21
                            if lvlf(cii)~=lvl(k)%bug fix JP3 10/4/21
                                B{kk}{rr,rrr+cii-1}=[];%bug fix JP3 10/4/21
                                y{kk}{rrr+cii-1}=[];%bug fix JP3 10/4/21
                            end%bug fix JP3 10/4/21
                        end%bug fix JP3 10/4/21
                        rrr=rrr+stoi-1;%bug fix JP3 10/4/21
                    end%bug fix JP3 10/4/21
                end%bug fix JP3 10/4/21
				if ~isempty(Aleak)
					B{kk}{rr,rrr}=[Aleak '*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
					y{kk}{rrr}=yAleak;
					rrr=rrr+1;
				end
% 				for cii=1:stoi
% 					if ~isempty(Bleak{cii})
% 						B{kk}{rr,rrr}=[Bleak{cii} '*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
% 						y{kk}{rrr}=[yleak{cii}];
% 						rrr=rrr+1;
% 					end
% 				end
			end
			for f=reverse{r}
				stoi=abs(minS(emu{r}.meti,f));
                stoichiometry=stoi; % stoichiometry is the absolute value of stoichiometry %bug fix JP 8/6/21
                rrr=rrr+1;
				Aleak=[];
				yAleak=[];
				if stoi==1
					A{kk}{rr,rr}=[A{kk}{rr,rr} '-f(' num2str(f) ')'];
				elseif stoi>1
					A{kk}{rr,rr}=[A{kk}{rr,rr} '-' num2str(stoi) '*f(' num2str(f) ')'];
				end
				prevf=zeros(I(k),4);
				lll=1;
				for c=find(cell2mat(cellfun(@(x) sum(x.meti==to{r}),emu,'UniformOutput',false)))' %c is the all emu indices which current r(th) emu generates
					if emu{c}.input==1 || c>=I(k+1) || ~isfield(rxnss{f},'rproduct')
						continue
					end
					cc=c-I(k)+1;
% 					stoi=abs(S(emu{r}.meti,f));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					ro=r;
					if emu{r}.equiv==0
                        stoi=1; % number of equivalent EMUs
						requiv=r;
					else
						requiv=find(cell2mat(cellfun(@(x) x.equiv==emu{r}.equiv,emu,'UniformOutput',false)))';
						stoi=length(requiv);
					end
					for r=requiv
						if length(rxnss{f}.rproduct)>1
							lp=find(cellfun(@(x) strcmp(x.id,emu{c}.id),rxnss{f}.rproduct));
						else
							lp=find(strcmp(rxnss{f}.rproduct.id,emu{c}.id));
						end
						if length(rxnss{f}.reduct)>1
							le=find(cellfun(@(x) strcmp(x.id,emu{r}.id),rxnss{f}.reduct));
						else
							le=find(strcmp(rxnss{f}.reduct.id,emu{r}.id));
	% 						le=strfind(rxnss{f}.reduct.id,emu{c}.id);
						end
						for llp=lp
							for lle=le
                                stoi2=0;%bug fix JP4 10/5/21
								if length(rxnss{f}.reduct)>1 && length(rxnss{f}.rproduct)>1
									if ismember(rxnss{f}.rproduct{llp}.cfg(emu{c}.cfg),rxnss{f}.reduct{lle}.cfg(emu{r}.cfg))
										if cc<1
											prevf(lll,1)=c;
											prevf(lll,2)=emu{c}.meti;
											prevf(lll,3)=sum(ismember(rxnss{f}.rproduct{llp}.cfg(emu{c}.cfg),rxnss{f}.reduct{lle}.cfg(emu{r}.cfg)));
											prevf(lll,4)=r;
											lll=lll+1;
                                        else
                                            if stoichiometry>1 && length(rxnss{f}.reduct{lle}.cfg)>=length(rxnss{f}.rproduct{llp}.cfg(emu{c}.cfg)) && sum(cellfun(@(x) strcmp(x.id,rxnss{f}.rproduct{llp}.id),rxnss{f}.rproduct))==stoichiometry && sum(cellfun(@(x) strcmp(x.id,rxnss{f}.reduct{lle}.id),rxnss{f}.reduct))==stoichiometry%bug fix JP4 10/5/21
                                                stoi2=stoichiometry;%bug fix JP4 10/5/21
                                                rxn(f)=stoi2;%bug fix JP5 10/5/21
%                                             end%bug fix JP4 10/5/21 %bug fix JP6 10/11/21
                                            elseif stoichiometry>1 && length(rxnss{f}.reduct{lle}.cfg)==length(rxnss{f}.rproduct{llp}.cfg) && sum(cellfun(@(x) strcmp(x.id,rxnss{f}.reduct{lle}.id),rxnss{f}.reduct))==stoichiometry%bug fix JP6 10/11/21
                                                stoi2=stoichiometry;%bug fix JP6 10/11/21
                                            end%bug fix JP6 10/11/21
											if stoi==1
                                                if stoichiometry==1 || stoi2>0 || (stoichiometry>1 && length(rxnss{f}.reduct{lle}.cfg)>length(rxnss{f}.rproduct{llp}.cfg(emu{c}.cfg))) || (stoichiometry>1 && length(rxnss{f}.reduct{lle}.cfg(emu{r}.cfg))<=length(rxnss{f}.rproduct{llp}.cfg)/stoichiometry) %bug fix JP 8/6/21 %bug fix JP4 10/5/21 %bug fix JP8 10/19/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')'];%bug fix JP 8/6/21
                                                else%bug fix JP 8/6/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+' num2str(stoichiometry) '*b(' num2str(f) ')'];%bug fix JP 8/6/21
                                                end%bug fix JP 8/6/21
											elseif stoi>1
                                                if stoichiometry==1 || stoi2>0 || (stoichiometry>1 && length(rxnss{f}.reduct{lle}.cfg)>length(rxnss{f}.rproduct{llp}.cfg(emu{c}.cfg))) || (stoichiometry>1 && length(rxnss{f}.reduct{lle}.cfg(emu{r}.cfg))<=length(rxnss{f}.rproduct{llp}.cfg)/stoichiometry) %bug fix JP 8/6/21 %bug fix JP4 10/5/21 %bug fix JP8 10/19/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')/' num2str(stoi)];%bug fix JP 8/6/21
                                                else%bug fix JP 8/6/21
                                                    A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')/' num2str(stoi/stoichiometry)];%bug fix JP 8/6/21
                                                end%bug fix JP 8/6/21
											end
										end
									end
								elseif length(rxnss{f}.reduct)>1 && length(rxnss{f}.rproduct)==1
									if [ismember(rxnss{f}.rproduct.cfg(emu{c}.cfg),rxnss{f}.reduct{lle}.cfg(emu{r}.cfg)) cc>0]
										if stoi==1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')'];
										elseif stoi>1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')/' num2str(stoi)];
										end
									end
								elseif length(rxnss{f}.reduct)==1 && length(rxnss{f}.rproduct)>1
									if ismember(rxnss{f}.rproduct{llp}.cfg(emu{c}.cfg),rxnss{f}.reduct.cfg(emu{r}.cfg))
										if cc<1
											prevf(lll,1)=c;
											prevf(lll,2)=emu{c}.meti;
											prevf(lll,3)=sum(ismember(rxnss{f}.rproduct{llp}.cfg(emu{c}.cfg),rxnss{f}.reduct.cfg(emu{r}.cfg)));
											prevf(lll,4)=r;
											lll=lll+1;
										else
											if stoi==1
												A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')'];
											elseif stoi>1
												A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')/' num2str(stoi)];
											end
										end
									end
								else
									if [ismember(rxnss{f}.rproduct.cfg(emu{c}.cfg),rxnss{f}.reduct.cfg(emu{r}.cfg)) cc>0]
										if stoi==1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')'];
										elseif stoi>1
											A{kk}{rr,cc}=[A{kk}{rr,cc} '+b(' num2str(f) ')/' num2str(stoi)];
										end
									end
								end
							end
						end
						if cc>0 && ~isempty(A{kk}{rr,cc})
							if A{kk}{rr,cc}(1)=='+'
								A{kk}{rr,cc}(1)=[];
							end
							if leak(f)>0
								leakcfg=find(ismember(rxnss{f}.rproduct.cfg,strrep(rxnss{f}.rproduct.cfg(emu{c}.cfg),rxnss{f}.leak{1},'')));
								if length(leakcfg)<lvl(k) && ~isempty(leakcfg)
% 									Aleak=A{kk}{rr,cc};
									Aleak=['-b(' num2str(f) ')'];
									A{kk}{rr,cc}=[A{kk}{rr,cc} '*(1-f(' num2str(length(leak)+sum(leak(1:f))) '))'];
									leakemu=cellfun(@(x) x.meti==emu{c}.meti & isequal(x.cfg,leakcfg),emu);
									yAleak=[',x' num2str(emu{leakemu}.network) '(' num2str(emu{leakemu}.index) '),x' num2str(emu{hemu}.network) '(' num2str(emu{hemu}.index) ')'];
								elseif isempty(leakcfg)
									hh=find(hemu)-I(k)+1;
									A{kk}{rr,cc}=[A{kk}{rr,cc} '*(1-f(' num2str(length(leak)+sum(leak(1:f))) '))'];
									A{kk}{rr,hh}=[A{kk}{rr,hh} '+b(' num2str(f) ')*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
									A{kk}{hh,cc}=[A{kk}{hh,cc} '+b(' num2str(f) ')*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
									A{kk}{hh,hh}=[A{kk}{hh,hh} '-b(' num2str(f) ')*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
									if A{kk}{rr,hh}(1)=='+'
										A{kk}{rr,hh}(1)=[];
									end
									if A{kk}{hh,cc}(1)=='+'
										A{kk}{hh,cc}(1)=[];
									end
% 									yAleak=[',x' num2str(emu{hemu}.network) '(' num2str(emu{hemu}.index) ')'];
								end
							end
						end
					end
					r=ro;
				end
				if emu{r}.equiv==0
					stoi=1; % number of equivalent EMUs
					requiv=r;
				else
					requiv=find(cell2mat(cellfun(@(x) x.equiv==emu{r}.equiv,emu,'UniformOutput',false)))';
					stoi=length(requiv);
				end
				lvlf=zeros(stoi,1);
% 				Bleak=cell(stoi,1);
% 				yleak=cell(stoi,1);
				if sum(sum(prevf)) && isempty(cell2mat(regexp(A{kk}(rr,:),['b\(' num2str(f) '\)(?=$|[^/])'])))
					[uprevf,~,ic]=unique(prevf(:,[2 4]),'stable','rows');
                    [~,ia]=unique(prevf,'rows');                            %bug fix JP2 10/1/21
                    iia=find(~ismember(1:max(ia),ia));                      %bug fix JP2 10/1/21
                    rxnfstep=1;%bug fix JP5 10/5/21
                    if rxn(f)~=1 && ~isempty(cell2mat(regexp(A{kk}(rr,:),[num2str(rxn(f)) '\*f\(' num2str(f) '\)(?=$|[^/])'])))%bug fix JP5 10/5/21
                        rxnfstep=rxn(f);%bug fix JP5 10/5/21
                    end%bug fix JP5 10/5/21
                    if rxnfstep==1%bug fix JP5 10/5/21
                    for iiia=iia%bug fix JP2 10/1/21
                        if prevf(iiia,3)>=lvl(k)-prevf(ia(end),3)%bug fix JP2 10/1/21
                            rc=find(requiv==prevf(iiia,4));%bug fix JP2 10/1/21
                            if prevf(iiia,3)<max(prevf(:,3)) %bug fix JP2 10/1/21
                                y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(iiia)}.network) '(' num2str(emu{prevf(iiia)}.index) ')'];%bug fix JP2 10/1/21
                                lvlf(rc)=lvlf(rc)+prevf(iiia,3);%bug fix JP2 10/1/21
                            end%bug fix JP2 10/1/21
                        end%bug fix JP2 10/1/21
                    end%bug fix JP2 10/1/21
                    end%bug fix JP4 10/5/21
					for iic=1:max(ic)
						if uprevf(iic,:)~=0
							mprevf=max(prevf(ic==iic,3));
							ci=find(ismember(prevf(:,2:4),[uprevf(iic,1) mprevf uprevf(iic,2)],'rows'));
							rc=find(requiv==uprevf(iic,2));
							for cii=1:rxnfstep:length(ci)%bug fix JP5 10/5/21
								if cii>1 && maxS(emu{prevf(ci(cii))}.meti,f)>1
									y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];
									lvlf(rc)=lvlf(rc)+mprevf;
								elseif cii==1
									y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];
									lvlf(rc)=lvlf(rc)+mprevf;
								end
							end
							if stoi==1
								B{kk}{rr,rrr+rc-1}=['-b(' num2str(f) ')'];
                                if rxnfstep>1%bug fix JP5 10/5/21
                                    B{kk}{rr,rrr+rc-1}=['-b(' num2str(f) ')*' num2str(rxn(f))];%bug fix JP5 10/5/21
                                end%bug fix JP4 10/5/21
							elseif stoi>1
								B{kk}{rr,rrr+rc-1}=['-b(' num2str(f) ')/' num2str(stoi)];
                                if rxnfstep>1%bug fix JP5 10/5/21
                                    B{kk}{rr,rrr+rc-1}=['-b(' num2str(f) ')/' num2str(stoi) '*' num2str(rxn(f))];%bug fix JP5 10/5/21
                                end%bug fix JP4 10/5/21
							end
% 							if leak(f)>0
% 								Bleak{rc}=B{kk}{rr,rrr+rc-1};
% 								B{kk}{rr,rrr+rc-1}=[B{kk}{rr,rrr+rc-1} '*(1-f(' num2str(length(leak)+sum(leak(1:f))) '))'];
% 							end
						end
                    end
                    
					for cii=1:stoi
						if lvlf(cii)~=lvl(k)
% 							Bleak{cii}=[];
% 							yleak{cii}=[];
							B{kk}{rr,rrr+cii-1}=[];
							y{kk}{rrr+cii-1}=[];
						end
					end
					rrr=rrr+stoi-1;
                end
                for cc=1:llvl%bug fix JP3 10/4/21
                    if sum(sum(prevf)) && ~isempty(cell2mat(regexp(A{kk}(rr,cc),['f\(' num2str(f) '\)(?=$|[^/])']))) && ~isempty(cell2mat(regexp(A{kk}(rr,cc),['b\(' num2str(f) '\)(?=$|[^/])'])))%bug fix JP3 10/4/21
                        [uprevf,~,ic]=unique(prevf(:,[2 4]),'stable','rows');%bug fix JP3 10/4/21
                        [~,ia]=unique(prevf,'rows');                            %bug fix JP3 10/4/21
                        iia=find(~ismember(1:max(ia),ia));                      %bug fix JP3 10/4/21
                        for iiia=iia%bug fix JP3 10/4/21
                            if prevf(iiia,3)>=lvl(k)-prevf(ia(end),3)%bug fix JP3 10/4/21
                                rc=find(requiv==prevf(iiia,4));%bug fix JP3 10/4/21
                                if prevf(iiia,3)<max(prevf(:,3))%bug fix JP3 10/4/21
                                    y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(iiia)}.network) '(' num2str(emu{prevf(iiia)}.index) ')'];%bug fix JP3 10/4/21
                                    lvlf(rc)=lvlf(rc)+prevf(iiia,3);%bug fix JP3 10/4/21
                                end%bug fix JP3 10/4/21
                            end%bug fix JP3 10/4/21
                        end%bug fix JP3 10/4/21
                        for iic=1:max(ic)%bug fix JP3 10/4/21
                            if uprevf(iic,:)~=0%bug fix JP3 10/4/21
                                mprevf=max(prevf(ic==iic,3));%bug fix JP3 10/4/21
                                ci=find(ismember(prevf(:,2:4),[uprevf(iic,1) mprevf uprevf(iic,2)],'rows'));%bug fix JP3 10/4/21
                                rc=find(requiv==uprevf(iic,2));%bug fix JP3 10/4/21
                                for cii=1:length(ci)%bug fix JP3 10/4/21
                                    if cii>1 && maxS(emu{prevf(ci(cii))}.meti,f)>1%bug fix JP3 10/4/21
                                        y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];%bug fix JP3 10/4/21
                                        lvlf(rc)=lvlf(rc)+mprevf;%bug fix JP3 10/4/21
                                    elseif cii==1%bug fix JP3 10/4/21
                                        y{kk}{rrr+rc-1}=[y{kk}{rrr+rc-1} ',x' num2str(emu{prevf(ci(cii))}.network) '(' num2str(emu{prevf(ci(cii))}.index) ')'];%bug fix JP3 10/4/21
                                        lvlf(rc)=lvlf(rc)+mprevf;%bug fix JP3 10/4/21
                                    end%bug fix JP3 10/4/21
                                end%bug fix JP3 10/4/21
                                if stoi==1%bug fix JP3 10/4/21
                                    B{kk}{rr,rrr+rc-1}=['-b(' num2str(f) ')'];%bug fix JP3 10/4/21
                                elseif stoi>1%bug fix JP3 10/4/21
                                    B{kk}{rr,rrr+rc-1}=['-b(' num2str(f) ')/' num2str(stoi)];%bug fix JP3 10/4/21
                                end%bug fix JP3 10/4/21
                            end%bug fix JP3 10/4/21
                        end%bug fix JP3 10/4/21
                        for cii=1:stoi%bug fix JP3 10/4/21
                            if lvlf(cii)~=lvl(k)%bug fix JP3 10/4/21
                                B{kk}{rr,rrr+cii-1}=[];%bug fix JP3 10/4/21
                                y{kk}{rrr+cii-1}=[];%bug fix JP3 10/4/21
                            end%bug fix JP3 10/4/21
                        end%bug fix JP3 10/4/21
                        rrr=rrr+stoi-1;%bug fix JP3 10/4/21
                    end%bug fix JP3 10/4/21
                end%bug fix JP3 10/4/21
				if ~isempty(Aleak)
					B{kk}{rr,rrr}=[Aleak '*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
					y{kk}{rrr}=yAleak;
					rrr=rrr+1;
				end
% 				for cii=1:stoi
% 					if ~isempty(Bleak{cii})
% 						B{kk}{rr,rrr}=[Bleak{cii} '*f(' num2str(length(leak)+sum(leak(1:f))) ')'];
% 						y{kk}{rrr}=[];
% 						rrr=rrr+1;
% 					end
% 				end
			end
		else
			rrr=rrr+1;
			B{kk}{cellfun(@(x) x.meti==to{r} && isequal(rxnss{reverse{r}}.rproduct.cfg(x.cfg),rxnss{reverse{r}}.reduct.cfg(emu{r}.cfg)),emu(I(k):I(k+1)-1)),rrr}=['-f(' num2str(reverse{r}) ')'];
			y{kk}{rrr}=['in.' emu{r}.id num2str(emu{r}.cfg,'_%d')];
%             if uexptype(mexpi)==1 || uexptype(mexpi)==3 || uexptype(mexpi)==7 || uexptype(mexpi)==15
            yy{yyy}=['in.' emu{r}.id num2str(emu{r}.cfg,'_%d') '=emufy(inp.' emu{r}.id '__,[' num2str(emu{r}.cfg,' %d') ']);'];
%             elseif uexptype(mexpi)==2
%             end
			yy{yyy}(strfind(yy{yyy},'  '))='';
			yy{yyy}(strfind(yy{yyy},' _'))='';
			yyy=yyy+1;
		end
	end
% reorganize A B and y; assign network and index
	newnet=1;
	maxnet=max(cell2mat(cellfun(@(x) x.network,emu,'UniformOutput',false)));
	emuind=1:size(A{kk},1);
	while sum(sum(~cellfun(@isempty,A{kk})))
		Alvl=~cellfun(@isempty,A{kk}(:,1));
		for i=1:size(A{kk},2)
			if sum(Alvl)==0
				Alvl=cellfun(@(x) ~isempty(x),A{kk}(:,i));
			end
			for j=1:size(A{kk},2)
				if sum(Alvl & cellfun(@(x) ~isempty(x),A{kk}(:,j)))
					Alvl=Alvl | cellfun(@(x) ~isempty(x),A{kk}(:,j));
				end
			end
		end
		A{kk+newnet}=A{kk}(:,Alvl);
		A{kk+newnet}=A{kk+newnet}(Alvl,:);
		A{kk}(:,Alvl)=[];
		A{kk}(Alvl,:)=[];
		emuind_new=emuind(Alvl)+I(k)-1;
		emuind(Alvl)=[];
		B{kk+newnet}=B{kk}(Alvl,:);
		B{kk}(Alvl,:)=[];
		y{kk+newnet}=y{kk};
		if size(B{kk+newnet},1)==1
			removerhs=logical(cellfun(@isempty,B{kk+newnet}));
		else
			removerhs=~logical(sum(cellfun(@(x) ~isempty(x),B{kk+newnet})));
		end
		B{kk+newnet}(:,removerhs)=[];
		y{kk+newnet}(removerhs)=[];
		newind=1;
		for newl=emuind_new
			emu{newl}.network=maxnet+newnet;
			emu{newl}.index=newind;
			newind=newind+1;
		end
		newnet=newnet+1;
	end
end

% clean up empty cells
yy=yy(~cellfun(@isempty,yy));
for i=length(lvl):-1:1
	for j=cl:-1:1
		if sum(sum(cellfun(@(x) ~isempty(x),A{cl*(i-1)+j})))==0
			A{cl*(i-1)+j}=[];
			B{cl*(i-1)+j}=[];
			y{cl*(i-1)+j}=[];
		end
	end
end

mexp{mexpi}.A=A(~cellfun(@isempty,A));
mexp{mexpi}.B=B(~cellfun(@isempty,B));
mexp{mexpi}.y=y(~cellfun(@isempty,y));
mexp{mexpi}.yy=yy;
mexp{mexpi}.emu=emu;

%% Print function to simulate EMU from free fluxes and input labeling
if isempty(strfind(optfilename,'.'))
	fid=fopen([optfilename '_' mexp{mexpi}.exptype '.m'],'w');
else
	fid=fopen([optfilename(1:strfind(optfilename,'.')-1) '_' mexp{mexpi}.exptype '_' optfilename(strfind(optfilename,'.')+1:end)],'w'); %%%%%%%%%%%%%%%%%%
	optfilename=optfilename(1:strfind(optfilename,'.')-1);
end

% % wrapper function to simulate multiple labeling experiments

currexp=find(exptypei'==mexpi);
for expi=currexp
	fprintf(fid_,'\npersistent in%s\nif isempty(in%s)\n',exp{expi},exp{expi});
	for i=1:length(yy)
		fprintf(fid_,'\t%s\n',[yy{i}(1:2) exp{expi} yy{i}(3:strfind(yy{i},'__')+1) exp{expi} yy{i}(strfind(yy{i},'__')+2:end)]);
	end
	fprintf(fid_,'end\n');
	fprintf(fid_,'mol_=%s_%s(free_net,free_xch,in%s);\n',optfilename,mexp{mexpi}.exptype,exp{expi});

	for i=1:length(mexp{mexpi}.metfrag)
		emuind=cellfun(@(x) x.meti==mexp{mexpi}.metfrag{i}.meti & isequal(x.cfg,mexp{mexpi}.metfrag{i}.cfg),emu);
		tempatoms=0;
		tempcfg=[];
		for j=find(de2bi(uexptype(mexpi),length(atoms)))
			tempatoms=tempatoms+metss{emu{emuind}.meti}.(atoms{j});
			if j==1
				tempcfg=1:tempatoms;
			else
				otheratoms=0;
				for k=1:j-1
					otheratoms=otheratoms+metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j-k});
				end
				tempcfg=[tempcfg otheratoms+(1:metss{emu{emuind}.meti}.(atoms{j}))];
			end;
		end
		if sum(emuind) && emu{emuind}.atoms==metss{emu{emuind}.meti}.atoms
			fprintf(fid_,'mol.%s__%s=mol_.%s;\n',emu{emuind}.id,exp{expi},emu{emuind}.id);
		elseif sum(emuind) && emu{emuind}.atoms<metss{emu{emuind}.meti}.atoms
			if isequal(emu{emuind}.cfg,tempcfg)
				fprintf(fid_,'mol.%s__%s=mol_.%s;\n',emu{emuind}.id,exp{expi},emu{emuind}.id);
			else
				fprintf(fid_,'mol.%s__%s',emu{emuind}.id,exp{expi});
				emucfg=num2str(emu{emuind}.cfg,'_%d');
				emucfg(strfind(emucfg,'  '))='';
				emucfg(strfind(emucfg,' _'))='';
				fprintf(fid_,'%s=mol_.%s%s;\n',emucfg,emu{emuind}.id,emucfg);
			end
		end
	end
end

% % universal simulation function with different input labeling
fprintf(fid,'function mol=%s_%s(free_net,free_xch,in)\n\n',optfilename,mexp{mexpi}.exptype);

fprintf(fid,'kernel_net=[');
for i=1:size(m.kernel_net,1)-1
	fprintf(fid,' %g',m.kernel_net(i,:));
	fprintf(fid,';\n');
end
fprintf(fid,' %g',m.kernel_net(i+1,:));
fprintf(fid,'];\n\n');

fprintf(fid,'kernel_xch=[');
for i=1:size(m.kernel_xch,1)-1
	fprintf(fid,' %g',m.kernel_xch(i,:));
	fprintf(fid,';\n');
end
fprintf(fid,' %g',m.kernel_xch(i+1,:));
fprintf(fid,'];\n\n');

fprintf(fid,'net=sparse(kernel_net)*free_net;\nxch=sparse(kernel_xch)*free_xch;\n\nf=xch+max(0,net);\nb=xch-min(0,net);\n');

for n=1:length(mexp{mexpi}.A)
	la=size(mexp{mexpi}.A{n},1);
	ly=length(mexp{mexpi}.y{n});
	
	fprintf(fid,'\n%%%%\nA%d=zeros(%d);\n',n,la);
	for i=1:la
		for j=1:la
			if ~isempty(mexp{mexpi}.A{n}{i,j})
				fprintf(fid,'A%d(%d,%d)=%s;\n',n,i,j,mexp{mexpi}.A{n}{i,j});
			end
		end
	end
	
	fprintf(fid,'\nB%d=zeros(%d,%d);\n',n,la,ly);
	for i=1:la
		for j=1:ly
			if ~isempty(mexp{mexpi}.B{n}{i,j})
				fprintf(fid,'B%d(%d,%d)=%s;\n',n,i,j,mexp{mexpi}.B{n}{i,j});
			end
		end
	end
	
	fprintf(fid,'\ny%d=[',n);
	for j=1:ly
		mexp{mexpi}.y{n}{j}(strfind(mexp{mexpi}.y{n}{j},'  '))='';
		mexp{mexpi}.y{n}{j}(strfind(mexp{mexpi}.y{n}{j},' _'))='';
		commaloc=strfind(mexp{mexpi}.y{n}{j},',');
		for i=length(commaloc)-1:-1:1
			mexp{mexpi}.y{n}{j}=[mexp{mexpi}.y{n}{j}(commaloc(1):commaloc(i)) 'cauchy(' mexp{mexpi}.y{n}{j}(commaloc(i)+1:end) ')'];
		end
		if mexp{mexpi}.y{n}{j}(1)==','
			mexp{mexpi}.y{n}{j}(1)=[];
		end
		commaloc=regexp(mexp{mexpi}.y{n}{j},'\d)');
		for i=length(commaloc):-1:1
			mexp{mexpi}.y{n}{j}=[mexp{mexpi}.y{n}{j}(1:commaloc(i)) ',:' mexp{mexpi}.y{n}{j}(commaloc(i)+1:end)];
		end
		if j<ly
			fprintf(fid,' %s;',mexp{mexpi}.y{n}{j});
		else
			fprintf(fid,' %s];',mexp{mexpi}.y{n}{j});
		end
	end
	
	if la>10
		fprintf(fid,'\nx%d=sparse(A%d)\\(B%d*sparse(y%d));\n',n,n,n,n);
	else
		fprintf(fid,'\nx%d=A%d\\(B%d*sparse(y%d));\n',n,n,n,n);
	end
end

fprintf(fid,'\n%%%%');
for i=1:length(mexp{mexpi}.metfrag)
	emuind=cellfun(@(x) x.meti==mexp{mexpi}.metfrag{i}.meti & isequal(x.cfg,mexp{mexpi}.metfrag{i}.cfg),emu);
	tempatoms=0;
	tempcfg=[];
	for j=find(de2bi(uexptype(mexpi),length(atoms)))
		tempatoms=tempatoms+metss{emu{emuind}.meti}.(atoms{j});
		if j==1
			tempcfg=1:tempatoms;
		else
			otheratoms=0;
			for k=1:j-1
				otheratoms=otheratoms+metss{strcmp(mets,metfrag_{i}.id)}.(atoms{j-k});
			end
			tempcfg=[tempcfg otheratoms+(1:metss{emu{emuind}.meti}.(atoms{j}))];
		end
	end
	if sum(emuind) && emu{emuind}.atoms==metss{emu{emuind}.meti}.atoms
		fprintf(fid,'\nmol.%s=x%d(%d,:);',emu{emuind}.id,emu{emuind}.network,emu{emuind}.index);
	elseif sum(emuind) && emu{emuind}.atoms<metss{emu{emuind}.meti}.atoms
		if isequal(emu{emuind}.cfg,tempcfg)
			fprintf(fid,'\nmol.%s=x%d(%d,:);',emu{emuind}.id,emu{emuind}.network,emu{emuind}.index);
		else
			fprintf(fid,'\nmol.%s',emu{emuind}.id);
			emucfg=num2str(emu{emuind}.cfg,'_%d');
			emucfg(strfind(emucfg,'  '))='';
			emucfg(strfind(emucfg,' _'))='';
			fprintf(fid,'%s=x%d(%d,:);',emucfg,emu{emuind}.network,emu{emuind}.index);
		end
	end
end
fclose(fid);
end
fclose(fid_);
% catch
% delete(h)
% warning('ERROR')
% end
