function  createStatsmap(vertex_struct,tempmgh,sub_num_threshold,distribution_descriptors, output_dir,activation_threshold)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = sub_num_threshold
    disp(['\nThis is subject threshold ' num2str(i)]);
    for m = 1: length(distribution_descriptors)
        metric= distribution_descriptors{m};
        funcmap= getFuncmap(vertex_struct, i, metric);
        % create the mgh
        file_store_dir= creatFuncmgh(funcmap,tempmgh,output_dir, i, metric, activation_threshold);
        disp(['\nThe output folder is : ' file_store_dir])
    end
end

end