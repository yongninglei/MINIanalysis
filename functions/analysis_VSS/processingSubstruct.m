function [subs_VOTfMRI,minmaxval] = processingSubstruct(VOTmaskedsubs,tempmgh,activation_threshold,processing_algo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fieldname= fieldnames(VOTmaskedsubs);
fieldname_f= fieldname(6:end);

% initialize struct for normalization 
% structure: 58x7x 163842
subs_VOTfMRI=VOTmaskedsubs;
delFields  = {'CT', 'CURV', 'qMRI_T1qMRI'};
subs_VOTfMRI = rmfield(subs_VOTfMRI, delFields);
NaNVector = NaN(size(tempmgh.vol'));

for ns = 1: length(VOTmaskedsubs)
    for i = 1: length(fieldname_f)
        fname = [fieldname_f{i}];
        subs_VOTfMRI(ns).(fname) = NaNVector;
    end
end

% get a intermediate struct stores the max value, this is for check
minmaxval=subs_VOTfMRI;
minmaxval=rmfield(minmaxval, {'VOT_block_RWvsCB','VOT_block_RWvsCS','VOT_block_RWvsFF','VOT_block_RWvsPS','VOT_block_RWvsSD','VOT_block_RWvsPW','VOT_block_RWvsNull'});
for sub_ind =1: length(subs_VOTfMRI)
   for i = 1: length(fieldname_f)
       fname = [fieldname_f{i}];
       fname_max = ['max_' fname];
       fname_min = ['min_' fname];
       minmaxval(sub_ind).(fname_max)    = max(VOTmaskedsubs(sub_ind).(fname));
       minmaxval(sub_ind).(fname_min)    = min(VOTmaskedsubs(sub_ind).(fname));
   end
end

% preprocessing the sub struct based on some criteiral, first all the value
% that are not valid will be transfered to nan, finally they will become 0
% for visulization
lex_filedname= {'VOT_block_RWvsCS','VOT_block_RWvsFF','VOT_block_RWvsPW'};
per_fieldname={'VOT_block_RWvsCB','VOT_block_RWvsPS','VOT_block_RWvsSD'};

for sub_ind =1: length(subs_VOTfMRI)
    for i = 1: length(fieldname_f)
        fname = [fieldname_f{i}];  
        array=VOTmaskedsubs(sub_ind).(fname);
        if strcmp(processing_algo, 'maxval_percentage_posonly')
            subs_VOTfMRI(sub_ind).(fname) = normalizeArray_percentage_posonly(array, activation_threshold);
        elseif strcmp(processing_algo, 'maxval_percentage')
            subs_VOTfMRI(sub_ind).(fname) = normalizeArray_percentage(array, activation_threshold);
        elseif strcmp(processing_algo, 'tavl_threshold')
            subs_VOTfMRI(sub_ind).(fname) = normalizeArray_tvalthreshold(array, activation_threshold);
        elseif strcmp(processing_algo, 'tavl_threshold_posonly')
            subs_VOTfMRI(sub_ind).(fname) = normalizeArray_tvalthreshold_posonly(array, activation_threshold);
        
        end
           
    end
    
    lex_array=[subs_VOTfMRI(sub_ind).('VOT_block_RWvsCS'),subs_VOTfMRI(sub_ind).('VOT_block_RWvsFF'),subs_VOTfMRI(sub_ind).('VOT_block_RWvsPW')];
    subs_VOTfMRI(sub_ind).('VOT_block_LEX')= mean(lex_array,2,'omitnan');
    per_array=[subs_VOTfMRI(sub_ind).('VOT_block_RWvsCB'),subs_VOTfMRI(sub_ind).('VOT_block_RWvsPS'),subs_VOTfMRI(sub_ind).('VOT_block_RWvsSD')];
    subs_VOTfMRI(sub_ind).('VOT_block_PER')= mean(per_array,2,'omitnan');
    
    new_fieldname= fieldnames(subs_VOTfMRI);
    filedname_func= new_fieldname(3:end);
    for ii = 1: length(filedname_func)
        fname_new=[filedname_func{ii}];
        % convert back all the nan value to 0
        subs_VOTfMRI(sub_ind).(fname_new)(isnan(subs_VOTfMRI(sub_ind).(fname_new))) =0;
    end
end

end