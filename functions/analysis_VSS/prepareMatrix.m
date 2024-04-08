function [matrix,col_name] = prepareMatrix(subs ,normalize)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
row_num= numel(subs);

% create a 0s matrix with the same row numbers and col numbers
matrix=zeros(row_num, numel(fieldnames(subs))-2);

filed_names= fieldnames(subs);
col_name=string(filed_names(3:end));

for i= 3:numel(filed_names)
    matrix(:,i-2)=[subs.(filed_names{i})];
end;
% normalize the matrix
if normalize
matrix=zscore(matrix);
end;

end