function batchMRIwritemgh(subidx,substruct,contrast,inputdir,outputdir, fsaverage_or_native)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
fname=sprintf('VOT_block_%s', contrast);
map=substruct(subidx).(fname);
sub=substruct(subidx).name;
file_store_dir= fullfile(outputdir,sub);

if ~exist(file_store_dir,'dir')
    mkdir(file_store_dir)
end

tempmgh= getTempmgh(inputdir,sub);
logP=tempmgh;
logP.vol= map;
if strcmp(fsaverage_or_native,'fsaverage')
    outputfile=[contrast '305.mgh'];
else
    outputfile=[contrast '.mgh'];
end

MRIwrite(logP,fullfile(file_store_dir, outputfile))
end