function [filtered_sub_struct, subject_index] = removeSubsFromStruct(subject_to_remove,struct_with_all_sub)
%This funtion is used to clear the struct based on an input cell array with
% 1xN size
%   Detailed explanation goes here


% get the index of all the subs that you need to remove, this index is used
% for LD.csv filter because that one is without subject name
removeInd= ismember({subject_to_remove.name}, subject_to_remove);

subject_index_to_remove= find(removeInd);
keep= true(size(subject_to_remove,1));
keep(subject_index_to_remove)=false;
ALLind=1:size(subject_to_remove,1);
subject_index= ALLind(keep);

filtered_sub_struct = subject_to_remove(subject_index);
end

