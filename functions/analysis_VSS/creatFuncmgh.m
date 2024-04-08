function file_store_dir=creatFuncmgh(funcmap, tempmgh,fMRI_threshold_dir, sub_num_threshold, measurement, activation_threshold)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

fieldname= fieldnames(funcmap);

file_store_dir= fullfile(fMRI_threshold_dir, ['activation_threshold_' num2str(activation_threshold)],...
    ['sub_num_threshold_',num2str(sub_num_threshold)], ['measurement_' measurement]);

if ~exist(file_store_dir,'dir')
    mkdir(file_store_dir)
end

for i = 1: length(fieldname)
    fname = [fieldname{i}];
    logP=tempmgh;
    logP.vol= [funcmap.(fname)];
    outputfile=[measurement '_' fname '_sub_' num2str(sub_num_threshold) '_activation_' num2str(activation_threshold) '_heatmap.mgh'];
    MRIwrite(logP,fullfile(file_store_dir, outputfile))
end



end