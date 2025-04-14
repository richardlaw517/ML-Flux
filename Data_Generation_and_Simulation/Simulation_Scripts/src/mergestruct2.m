function [singlestruct,singlestruct_uniform_w,maxw] = mergestruct2(structmatrix)
%% mergestruct merges struct vector into a single struct
% input: a vector of struct elements with fields that are 2D matrices.
% output: a single struct with elements that are concatenated matrices
% (with each metabolite having natural MID or uniform-length MID); and maxw
% the length of the largest metabolite MID vector.

fields=fieldnames(structmatrix);
maxw=1;
for i=1:numel(fields)
    w=size(structmatrix(1).(fields{i}),2);
    if w>maxw
        maxw=w;
    end
    for j=length(structmatrix):-1:1
        singlestruct.(fields{i})(:,(j-1)*w+1:j*w)=structmatrix(j).(fields{i});
    end
end
if nargout>1
    for i=1:numel(fields)
        w=size(structmatrix(1).(fields{i}),2);
        j=length(structmatrix);
        singlestruct_uniform_w.(fields{i})(:,(j-1)*maxw+1:(j-1)*maxw+w)=structmatrix(j).(fields{i});
        singlestruct_uniform_w.(fields{i})(:,(j-1)*maxw+w+1:j*maxw)=0;
        for j=length(structmatrix)-1:-1:1
            singlestruct_uniform_w.(fields{i})(:,(j-1)*maxw+1:(j-1)*maxw+w)=structmatrix(j).(fields{i});
        end
    end
end