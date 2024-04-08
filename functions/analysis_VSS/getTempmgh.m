function tempmgh = getTempmgh(datadir,sub)
%getTempmgh is used to get template mgh file to be used in MRIwrite
%   At the individual subejct level analysis, it needs to be one for each
%   subejct
    tempmgh = MRIread([datadir fsp sub fsp 'results' fsp 'AllNoFacesvsNull305.mgh']);
    tempmgh.vol = zeros(size(tempmgh.vol));
end