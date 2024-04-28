function [filtered_sub_struct, subject_index] = getSubStruct(subject_to_remove,subs,TESTind)
%This funciton is used to  get a subset of subejuects according to two
%things: the index of subjects that you want to have and the names of
%subjects you want to disgard
%   Detailed explanation goes here

% first get all the subjects in the subs struct
names= {subs.name};
% get the index of all the subs that you need to remove, this index is used
% for LD.csv filter because that one is without subject name
removeInd= ismember({subs.name}, subject_to_remove);

subject_index_to_remove= find(removeInd);
keep= true(size(TESTind));
keep(subject_index_to_remove)=false;
subject_index= TESTind(keep);

filtered_sub_struct = subs(subject_index);
end

