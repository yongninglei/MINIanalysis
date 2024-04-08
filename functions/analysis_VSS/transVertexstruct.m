function [vertex_struct,vertex_valuecount] = transVertexstruct(subs_VOTfMRI,tempmgh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fieldname= fieldnames(subs_VOTfMRI);
fieldname_f= fieldname(3:end);


%% initialize the new matrix
% with 163842 x7 x 58
clear vertex_struct
row_num=length(tempmgh.vol');
vertex_struct(row_num)= struct();
zerocell=zeros(length(subs_VOTfMRI),1);
for i = 1: length(fieldname_f)
    fname = [fieldname_f{i}];     
    [vertex_struct(:).(fname)] = deal(zerocell); 
end


for vertex_idx=1:length(vertex_struct)
    for i = 1: length(fieldname_f)
        fname = [fieldname_f{i}];
        subsarray=zeros(length(subs_VOTfMRI),1);
        for subid = 1: length(subs_VOTfMRI)
            subsarray(subid)=subs_VOTfMRI(subid).(fname)(vertex_idx);
            vertex_struct(vertex_idx).(fname)=subsarray;
        end
    end
end

% check how many non-zeros are there for all the vertex
vertex_valuecount=vertex_struct;

for vertex_idx=1:length(vertex_struct)
    for i = 1: length(fieldname_f)
        fname = [fieldname_f{i}];
        vertex_valuecount(vertex_idx).(fname)= nnz(vertex_struct(vertex_idx).(fname));
    end
end

end