function vertex_struct = getVertexStruct(index,subs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

field_names=fieldnames(subs);
vertex_struct=subs;

for i =1:numel(subs)
    for j=3:numel(field_names)
        field_name=field_names{j};
        vertex_struct(i).(field_name)=subs(i).(field_name)(index);
    end;
end;

end