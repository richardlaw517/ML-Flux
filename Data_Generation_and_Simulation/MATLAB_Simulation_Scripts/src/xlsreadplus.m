function [org_avg,covar,org_raw,avgvar] = xlsreadplus(xlsname,sheets,model)
%% Reads excel file containing labeling data and create a struct
% Inputs: xlsname (excel file name), sheets (cell array with names of the sheets), and model (structure with the model)
% Outputs: org_avg (average labeling patterns), covar (covariance matrix),
% org_raw (raw labeling patterns), and avgvar (average variance)

% NOTE: May need to add in error checking mechanism (e.g. correct # carbons)

% Check if model is included 
if nargin < 3 
    model = 0;
end

avgvar = 0;
num_sheets = length(sheets);
nmea = zeros(num_sheets+1,1);
covmat = cell(num_sheets,1);

for i = 1:num_sheets
    if ispc % Logical for if this is being run on Windows (pc)
        [ndata,metabolites]=xlsread(xlsname,sheets{i});
    elseif isunix % Logical for if this is being run on Linus or MacOS platform
        warning('off','MATLAB:xlsread:Mode')
        [ndata,metabolites]=xlsread(xlsname,sheets{i},'','basic');
    else
        disp('Operating system not recognized/supported (i.e.: not PC, Linus, or MacOS')
    end

    if isnan(sum(ndata,'all'))
        warning('Covariance matrix contains NaN')
    end

    ndata = ndata';
    metabolites = metabolites(:,1);
    I=find(~cellfun(@isempty,metabolites));
    I(end+1) = 1+size(ndata,2);

    % Removes any metabolite data that is not in the model
    if isstruct(model)
        for j = size(I,1)-1:-1:1
            metaboliteName = metabolites{I(j)};
            underscoreIndex = strfind(metaboliteName,'__');
            metaboliteName = metaboliteName(1:underscoreIndex-1);
            if ~any(strcmpi(model.mets,metaboliteName)) % If the metabolite is not in the model
                if j == size(I,1)-1
                    metabolites(I(j)) = [];
                else
                    metabolites(I(j):I(j+1)-1) = [];
                end
                ndata(:,I(j):I(j+1)-1) = [];
            end
        end
    end
    
    covmat{i}=nancov(ndata,'pairwise');
	if isscalar(covmat{i})
		covmat{i}=0;
	end
    avgvar=avgvar+sum(diag(covmat{i}));
    nmea(i+1)=nmea(i)+size(ndata,2);

    % Removes any non-existent data (i.e. NaN)
    I=find(~cellfun(@isempty,metabolites));
    I(end+1)=1+size(ndata,2);
    for j=1:length(I)-1
        org_raw.(metabolites{I(j)})=ndata(:,I(j):I(j+1)-1);
        nanrow=find(isnan(org_raw.(metabolites{I(j)})(:,1)));
        for k=1:length(nanrow)
            org_raw.(metabolites{I(j)})(nanrow(k),:)=[];
            nanrow=nanrow-1;
        end
        org_avg.(metabolites{I(j)})=mean(org_raw.(metabolites{I(j)}),1);
    end

    covar(nmea(i)+1:nmea(i+1),nmea(i)+1:nmea(i+1))=covmat{i};
end

% average variance if not using covariance matrix
avgvar=avgvar/nmea(end);
if isnan(avgvar)
    warning('If every sample has NaN, the row should be removed.')
end

end
