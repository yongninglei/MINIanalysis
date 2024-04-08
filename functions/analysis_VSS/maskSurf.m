function maskedsubs = maskSurf(TESTsubs,label,tempmgh)
% 
% This function is used to using a label to mask the surface
% 
% The input would be: 
%   
%   TESTsubs: a struct in which each cell is a 163842x1 array, col are the
%   measurement and row are the subject
%   
%   label: fs label fule
%   
%   tempmgh: a template mgh to be used as a reference to write
% 
% Output: 
%   
%   A struct with the same dimension but each cell is masked by the label
%
% Copyright Yongning Lei (BCBL t.lei@bcbl.eu)


%% Apply VOTC mask to all filed and then save
maskedsubs=TESTsubs;
fieldN=fieldnames(TESTsubs);
for i =1:numel(TESTsubs)
    for j=3:numel(fieldN)
        fieldname=fieldN{j};
        allvertex    = [TESTsubs(i).(fieldname)]';
        % this step is using the lable to mask out the part we don't want
        allvertex    = allvertex(:, label);
        masked  = repmat(tempmgh.vol, [size(allvertex,1),1]);
        masked(label)=allvertex;
        maskedsubs(i).(fieldname)=masked';
    end;
end;
end