function tempmgh = getTempmgh(datadir,sub)
%getTempmgh is used to get template mgh file to be used in MRIwrite
%   At the individual subejct level analysis, it needs to be one for each
%   subejct
fsp=filesep;    
tempmgh = MRIread([datadir fsp sub fsp 'results' fsp 'AllNoFacesvsNull.mgh']);
tempmgh.vol = zeros(size(tempmgh.vol));

end