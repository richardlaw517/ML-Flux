function singlestruct = mergestruct(structvector)
%% mergestruct merges struct vector into a single struct
% input: a vector of struct elements with fields that are vectors
% output: a single struct with elements that are matrices

fields=fieldnames(structvector);
for i=1:numel(fields)
    for j=length(structvector):-1:1
        singlestruct.(fields{i})(j,:)=structvector(j).(fields{i});
    end
end