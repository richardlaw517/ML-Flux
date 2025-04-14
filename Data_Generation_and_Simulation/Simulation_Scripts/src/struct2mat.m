function catmat = struct2mat(structmat,dim)
%% struct2mat concatenates matrices in struct into a single matrix
% input: a struct with fields that are matrices (or vectors) and dimension
% along which matrices are concatenated (vertically if dim==1; horizontally
% if dim=2). Default dim is 2
% output: a single concatenated matrix

if nargin==1 || isempty(dim)
    dim=2;
end

fields=fieldnames(structmat);
catmat=structmat.(fields{2});
for i=3:numel(fields)
    catmat=cat(dim,catmat,structmat.(fields{i}));
end