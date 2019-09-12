function [H,P,CI,STATS] = myqMRIAnalysisInROIs(trt, subs, kkvertex, tempmgh)
% myqMRIAnalysisInROIs('TEST', TESTsubs, tempmgh)
    % trt = 'TEST'
    % subs = TESTsubs
    % trt = 'RETEST'
    % subs = DAY2subs

T1    = [subs.qMRI_T1qMRI]';

% T1 - litVWFA
% T1a = T1(:,kkvertex.aVWFA);
% T1c = T1(:,kkvertex.cVWFA);
% T1p = T1(:,kkvertex.pVWFA);
% mT1a    = mean(T1a, 2, 'omitnan');
% mT1c    = mean(T1c, 2, 'omitnan');
% mT1p    = mean(T1p, 2, 'omitnan');
% [H,P,CI,STATS] = ttest(mT1a, mT1c, 'alpha',0.05,'dim',1,'tail','both');
% [H,P,CI,STATS] = ttest(mT1a, mT1p, 'alpha',0.05,'dim',1,'tail','both');
% [H,P,CI,STATS] = ttest(mT1c, mT1p, 'alpha',0.05,'dim',1,'tail','both');

% T1 - fMRI results
T1per4 = T1(:,kkvertex.([trt '_perVWFA4']));
T1lex4 = T1(:,kkvertex.([trt '_lexVWFA4']));
mT1per4 = mean(T1per4, 2, 'omitnan');
mT1lex4 = mean(T1lex4, 2, 'omitnan');

[H,P,CI,STATS] = ttest(mT1per4, mT1lex4, 'alpha',0.05,'dim',1,'tail','both');
% csvwrite(fullfile(MINIPath,'DATA','qMRI', [trt '_mT1per4_T1.csv']),mT1per4);
% csvwrite(fullfile(MINIPath,'DATA','qMRI', [trt '_mT1lex4_T1.csv']),mT1lex4);

% Create table in long format and add subject information to be more clear in R
datos = array2table([mT1per4;mT1lex4]);
datos.Properties.VariableNames = {'T1'};
tmpsubs = struct2cell(subs);
tmpsubs = [tmpsubs(1,:),tmpsubs(1,:)];
datos.SUBJECT = tmpsubs';
datos.VWFA    = [repmat({'pOTS'},[height(datos)/2, 1]) ; ...
                 repmat({'mOTS'},[height(datos)/2, 1])];
             
writetable(datos, ...
           fullfile(MINIPath,'DATA','qMRI', [trt '_mOTS_pOTS_T1.csv']), ...
           'FileType', 'text', ...
           'Delimiter', 'comma', ...
           'WriteVariableNames', true)
           



% Save the T1 values for fsaverage visualization for Figure 5
T1aacacppp = T1;
T1aacacppp(:) = 0;
T1aacacppp(:,kkvertex.aacacppp) = T1(:,kkvertex.aacacppp);
idatzi = tempmgh;
idatzi.vol = mean(T1aacacppp,1,'omitnan');
MRIwrite(idatzi, fullfile(MINIPath,'DATA','qMRI', [trt '_T1aacacppp.mgh']));  

end