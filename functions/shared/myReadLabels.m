function kkvertex = myReadLabels(labeldir)

fsp = filesep;

% For the probabilistic fMRI analysis
adjLOTS = myFSread_label('fsaverage',[labeldir fsp 'lh.adjLOTS2.label'], 1);
kkvertex.adjLOTS = adjLOTS(:,1) + 1;

% Rosenke Cytoarchitectonic maps
rosenkelabels = dir([labeldir fsp 'Rosenke_vcAtlas_FS' fsp '*.label']);
for nl = 1: length(rosenkelabels)
    tmplabel = myFSread_label('fsaverage',...
                              [labeldir fsp 'Rosenke_vcAtlas_FS' fsp ...
                               rosenkelabels(nl).name], 1);
    myname = strrep(rosenkelabels(nl).name,'.label','');
    myname = strrep(myname,'.','_');
    kkvertex.(myname) = tmplabel(:,1) + 1;
end

% aparc labels
aparcs = dir([labeldir fsp 'aparcLabels' fsp '*.label']);
for nl = 1: length(aparcs)
    tmplabel = myFSread_label('fsaverage',...
                              [labeldir fsp 'aparcLabels' fsp ...
                               aparcs(nl).name], 1);
    myname = strrep(aparcs(nl).name,'.label','');
    myname = strrep(myname,'.','_');
    kkvertex.(myname) = tmplabel(:,1) + 1;
end

% mOTS and pOTS
TEST_perVWFA4 = myFSread_label('fsaverage',[labeldir fsp 'TEST_VOT_b_Perceptual4.label'], 1);
TEST_lexVWFA4 = myFSread_label('fsaverage',[labeldir fsp 'TEST_VOT_b_Lexical4.label'], 1);
kkvertex.pOTS = TEST_perVWFA4(:,1) + 1;
kkvertex.mOTS = TEST_lexVWFA4(:,1) + 1;

% Read V1 & BRODMANNS
V1 = myFSread_label('fsaverage',[labeldir fsp 'lh.V1.thresh.label'], 1);
kkvertex.V1 = V1(:,1) + 1;
V2 = myFSread_label('fsaverage',[labeldir fsp 'lh.V2.thresh.label'], 1);
kkvertex.V2 = V2(:,1) + 1;
BA44 = myFSread_label('fsaverage',[labeldir fsp 'lh.BA44.thresh.label'], 1);
kkvertex.BA44 = BA44(:,1) + 1;
BA45 = myFSread_label('fsaverage',[labeldir fsp 'lh.BA45.thresh.label'], 1);
kkvertex.BA45 = BA45(:,1) + 1;

% Read the fROIS
Anova_Groups_48_61_14 = myFSread_label('fsaverage',[labeldir fsp 'fROI' fsp 'Anova_Groups_-48_-61_-14.label'], 1);
kkvertex.Anova_Groups_48_61_14 = Anova_Groups_48_61_14(:,1) + 1;
Left_InfOper_Masked_AllControl_FWE01_10k = myFSread_label('fsaverage',[labeldir fsp 'fROI' fsp 'Left_InfOper_Masked_All-Control_FWE01_10k.label'], 1);
kkvertex.Left_InfOper_Masked_AllControl_FWE01_10k = Left_InfOper_Masked_AllControl_FWE01_10k(:,1) + 1;






end