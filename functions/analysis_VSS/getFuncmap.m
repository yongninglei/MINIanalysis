function vertex_final = getFuncmap(vertex_struct,sub_num_threshold,measurement)

fieldname= fieldnames(vertex_struct);

vertex_final= vertex_struct;
for vertex_idx=1:length(vertex_struct)
    for i = 1: length(fieldname)
        fname = [fieldname{i}];
        num_of_sub= nnz(vertex_struct(vertex_idx).(fname));
        if num_of_sub<sub_num_threshold
            vertex_final(vertex_idx).(fname)=0 ;
        else
            array=vertex_struct(vertex_idx).(fname);
            non_zero_array= array(array~=0);

            if strcmp(measurement,'mean')
                vertex_final(vertex_idx).(fname)=mean(non_zero_array) ;
            elseif strcmp(measurement,'median')
                vertex_final(vertex_idx).(fname)=median(non_zero_array);
            elseif strcmp(measurement,'std')
                vertex_final(vertex_idx).(fname)=std(non_zero_array);
            end
        end
    end
end

end