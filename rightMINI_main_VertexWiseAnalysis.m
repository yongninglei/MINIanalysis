%% Set up folders and variables
% Script based on the exploratory analysis file ~/code/MATLAB/myGlmfitPerVertex.m
% It has been reduced to contain only the code needed to reproduce the PNAS
% paper results.
% TODO1 = Use tables instead of structs... 
clc; clear all;
%tbUse paper-rightMINI;
tbUse MINIanalysis
fsp = filesep;


% Folders
MINIanalysis= '/bcbl/home/home_n-z/tlei/toolboxes/MINIanalysis';
MINIdir = '/bcbl/home/public/Gari/MINI';
% where the data are stoted, point to a directory with a lot of subjects
DATAdir = [MINIdir fsp 'ANALYSIS'];

% fMRI result subject dir
fmriresultdir= '/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/block/analysis_block_acpc_lhITfusLatOcc/SUBJECTS'
ANALYSISdir = [MINIdir fsp 'ANALYSIS' fsp 'fs_glmfit'];
fsdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc';
fMRI_threshold_dir= [ANALYSISdir fsp 'fMRI_threshold'] ;
% remember to setenv!!! 
% those two are not needed in lmc02
%fsbin = '/Applications/freesurfer/bin';
%fshome = '/Applications/freesurfer'; 
% labeldir = '/Applications/freesurfer/subjects/fsaverage/label';  
% labeldiraparc = '/Applications/freesurfer/subjects/fsaverage/label/aparcLabels';
% This was very very bad, I lost all myLabels when I upgraded the laptop, because I did
% not backup the Applications folder, because I thought no data was there. 

% In the server there is no this problem, only locally.

labeldir      = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc/fsaverage/label';

% Options
% Select to read
smmm = '5';
sm = '.fhmw5'; % '.fhmw5': Read smoothed data, '' read non smoothed data
SHOW = 0;  %  1: Verbose, 0: Quiet


% Separate in subjects in groups
% Consider that: 
% ---- 'S029', 'S032' 'S063' : do not have LD behaviroal result
% ---- 'S067' has no DWI
% ---- 'S013', 'S018' : do not use for fMRI BLOCK (didn't pass QA)
% ---- 'S013', 'S072', 'S075': do not use for fMRI EVENT (didn't pass QA)
% ---- 'S004', 'S029', 'S032','S048', 'S056', 'S072', 'S086': missing qMRI data
% ---- 'S067', 'S097' is the Polish girl
ALLind  = [1:97];
ONCEind = [1:35];
DAY1ind = [36:65,96];
DAY2ind = [66:95,97];
TESTind = [ONCEind, DAY1ind];

% Lists used to iterate afterwards
fMRIareas = {'VOT'};
designs   = {'block'};
ContrastSinNull = {'RWvsCB','RWvsCS','RWvsFF','RWvsPS','RWvsPW','RWvsSD'};
% I want to add faces and create a new contrast group, to see how the
% cluster goes, should be funny
ContrastConNull = {'RWvsCB','RWvsCS','RWvsFF','RWvsPS','RWvsPW','RWvsSD','RWvsNull'};
Contrasts       = ContrastConNull;
 

% Read all subjects in the folder and delete the non interesting fields of the
% struct
cd('/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/block/analysis_block_acpc_lhITfusLatOcc/smoothed_fhmw5');
subs = dir('S*');
delFields  = {'date', 'bytes', 'isdir', 'datenum'};
subs = rmfield(subs, delFields);
% Read a reference mgh to be used afterwards
tempmgh = MRIread([fmriresultdir fsp 'S001' fsp 'results' fsp 'AllNoFacesvsNull305.mgh']);
tempmgh.vol = zeros(size(tempmgh.vol));

%% Initialize subject structures and read them 
% Initialize structure with NaNs, and fill with data if it exists
% Now I can read qMRI MTV, T1 functional as well as curv and CT with sm 5
subs = myInitStructures(subs, tempmgh, fMRIareas, designs, Contrasts);

% SUPPORT DATA CREATION (just once, then leave it commented)
% Smoothear surface data 
% mySmoothStructures(subs, smmm, DATAdir, SHOW,'lh')
% Create labels with info coming from R analyses
% myCreateLabelsfsaverage('TEST')
% myCreateLabelsfsaverage('RETEST')
%% for tiger testing, needs to go to the function and read the data by hand
% Read MRI data coming from the server analyses
% for tiger testub
subs = myReadStructures(subs, DATAdir, sm, tempmgh, SHOW,'lh');
% Read the behavioural data
LD = readtable(['/bcbl/home/public/Gari/MINI/DATA' fsp 'LD.csv']); 

% Read the labels created in this project, fs305 space, analysis specific
% When reading lable files, remember to add 1 to the index 
kkvertex = myReadLabels(DATAdir, labeldir, ANALYSISdir);

%% Clean data according to QA
% % Make NaN those subjects that didn't pass QA
% NaNVector = NaN(size(tempmgh.vol'));
% % ---- 'S067' has no DWI >>> subs(62).name, and it is a bad subject
% % ---- 'S097' is the second day, make it nan as well: subs(92).name
% subs(62).DWI_vOF             = NaNVector;
% subs(62).DWI_pARC            = NaNVector;
% subs(92).DWI_vOF             = NaNVector;
% subs(92).DWI_pARC            = NaNVector;
% % ---- 'S013', 'S072', 'S075': do not use for fMRI EVENT (didn't pass QA)
% nsrm = [13, 67, 70];
% for ns = nsrm
%     for area = fMRIareas; for design = {'event'}; for contrast = ContrastSinNull
%         fname = [area{:} '_' design{:} '_' contrast{:}];
%          subs(ns).(fname)    = NaNVector;
%     end; end; end
% end
% % ---- 'S013', 'S018' : do not use for fMRI VOT,IFG,PPC BLOCK (didn't pass QA)
% nsrm = [13, 18];
% for ns = nsrm
%     for area = fMRIareas; for design = {'block'}; for contrast = ContrastSinNull
%         fname = [area{:} '_' design{:} '_' contrast{:}];
%          subs(ns).(fname)    = NaNVector;
%     end; end; end
% end
% ONCEsubs = subs(ONCEind);
% DAY1subs = subs(DAY1ind);
% # Write down why this subjects have to be removed from the analysis
% # S072: qMRI missing, the pipeline did not end, usually registration problems
% # S086: qMRI missing, the pipeline did not end, usually registration problems
% # S097: el RETEST de lexical decision
% # S075:  no event related data
% TESTsubs = subs(TESTind);
% DAY2subs = subs(DAY2ind);
% ALLsubs  = subs;
% subject_index = ALLind;

% right now didn't touch the QA of qMRI and DWI so I tick out the
% participants

%trt= 'TEST';

subjects_to_remove={'S004', 'S013', 'S018', 'S029', 'S032', 'S048', 'S056', 'S067' };

[TESTsubs, subject_index] = getSubStruct(subjects_to_remove,subs,TESTind);

%% Apply VOTC mask to all filed and then save
% run them once, now we just load
VOTmaskedsubs_load=load([ANALYSISdir fsp 'VOTmaskedsubs_CURVCTqMRIT1qMRIMTV6fMRI.mat']).VOTmaskedsubs;
% VOTmaskedsubs=maskSurf(TESTsubs,kkvertex.('VOT'),tempmgh);
% and then save it 
% save([ANALYSISdir fsp 'VOTmaskedsubs_CURVCTqMRIT1qMRIMTV6fMRI.mat'],'VOTmaskedsubs')
%% maintain everything but seperate PER and LEX from the struct

LEXsubs=load([ANALYSISdir fsp 'VOTmaskedsubs_LEX_CURVCTqMRIT1qMRIMTV.mat']).LEXsubs;
PERsubs=load([ANALYSISdir fsp 'VOTmaskedsubs_PER_CURVCTqMRIT1qMRIMTV.mat']).PERsubs;

% PERsubs=VOTmaskedsubs;
% LEXsubs=VOTmaskedsubs;
% 
% LEXdelFields  = {'VOT_block_RWvsCB','VOT_block_RWvsPS','VOT_block_RWvsSD','VOT_block_RWvsNull' };
% PERdelFields  = {'VOT_block_RWvsCS','VOT_block_RWvsFF','VOT_block_RWvsPW','VOT_block_RWvsNull' };
% LEXsubs = rmfield(LEXsubs, LEXdelFields);
% PERsubs = rmfield(PERsubs, PERdelFields);
% 
% save([ANALYSISdir fsp 'VOTmaskedsubs_LEX_CURVCTqMRIT1qMRIMTV.mat'],'LEXsubs');
% save([ANALYSISdir fsp 'VOTmaskedsubs_PER_CURVCTqMRIT1qMRIMTV.mat'],'PERsubs');
%% get some vertex for testing, get 1 per and get 1 sub 

index=162119;
vertexStruct_LEX=getVertexStruct(index,LEXsubs);
vertexStruct_PER=getVertexStruct(index,PERsubs);

index=2489;
vertex_2489_Struct_LEX=getVertexStruct(index,LEXsubs);
vertex_2489_Struct_PER=getVertexStruct(index,PERsubs);
%save([ANALYSISdir fsp 'vertex_162119_LEX.mat'],'vertexStruct_LEX');
%save([ANALYSISdir fsp 'vertex7_PER.mat'],'vertexStruct_PER');
%% Do the pure fMRI plot and check the pattern
% load the data, list all the parameters that we need to use
fname='VOTmaskedsubs_CURVCTqMRIT1qMRIMTV6fMRI.mat';
VOTmaskedsubs=load(fullfile(ANALYSISdir,fname)).VOTmaskedsubs;
measurements={'mean', 'median'};
%threshs={0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85};
threshs={0, 0.3,0.5};
output_dir=fullfile(fMRI_threshold_dir, 'non_negative_val');
processing_algo='maxval_percentage';  % valid choise 'no_negative_maxval_percentage'; 'maxval_percentage' 'tavl_threshold';

for at = 1:length(threshs) 
    activation_threshold= threshs{at};
    disp(['\nThis is normalization threshold ' num2str(activation_threshold)]);
    % normalize subs struct based on every's min max T val and turn that into a
    % -1 to 1 statistical struct
    [subs_VOTfMRI,vertex_struct,vertex_valuecount,minmaxval] = getVertexstruct(VOTmaskedsubs,tempmgh,activation_threshold,processing_algo);
    % threshold at the group level
    %{
    for every vertex, I will do a check, if this vertex have more than 20
    subjects that are activate, then I will say this vertex is 
    %}
    for i = 1:40
        disp(['\nThis is subject threshold ' num2str(i)]);
        for m = 1: length(measurements)
            measurement= measurements{m};
            funcmap= getFuncmap(vertex_struct, i, measurement);
            % create the mgh
            file_store_dir= creatFuncmgh(funcmap,tempmgh,output_dir, i, measurement, activation_threshold);
            disp(['\nThe output folder is : ' file_store_dir])
        end 
    end
end

% generate std for all the voxel when seeing different stimuli
activation_threshold= 0;
output_dir=fullfile(fMRI_threshold_dir, 'non_negative_val');
processing_algo='maxval_percentage';  % valid choise 'maxval_percentage_posonly'; 'maxval_percentage' 'tavl_threshold';
disp(['This is normalization threshold ' num2str(activation_threshold)]);
[subs_VOTfMRI,vertex_struct,vertex_valuecount,~] = getVertexstruct(VOTmaskedsubs,tempmgh,activation_threshold,processing_algo);

for i = 8:40
    disp(['\nThis is subject threshold ' num2str(i)]);
    measurement= 'std';
    funcmap= getFuncmap(vertex_struct, i, measurement);
    % create the mgh
    file_store_dir= creatFuncmgh(funcmap,tempmgh,output_dir, i, measurement, activation_threshold);
    disp(['The output folder is : ' file_store_dir])
end
%% generate T VAL 1.65 threshold result
% threshold all perticipant using CV=1.65, this represnets the
% significantly valid voxel for each subject.

% load the data, list all the parameters that we need to use
fname='VOTmaskedsubs_CURVCTqMRIT1qMRIMTV6fMRI.mat';
VOTmaskedsubs=load(fullfile(ANALYSISdir,fname)).VOTmaskedsubs;
distribution_descriptors={'mean', 'median','std'};
output_dir=fullfile(fMRI_threshold_dir, 'tval_0_threshold_posonly');
processing_algo='tavl_threshold_posonly';
tval_threshold=0; % before it's 1.65
sub_num_threshold=1:29;

[subs_VOTfMRI,minmaxval] = processingSubstruct(VOTmaskedsubs,tempmgh,activation_threshold,processing_algo);
[vertex_struct,vertex_valuecount] = transVertexstruct(subs_VOTfMRI,tempmgh);

createStatsmap(vertex_struct,tempmgh,sub_num_threshold,distribution_descriptors, output_dir,tval_threshold)


%% aggregate LEX and PER
fname='VOTmaskedsubs_CURVCTqMRIT1qMRIMTV6fMRI.mat';
VOTmaskedsubs=load(fullfile(ANALYSISdir,fname)).VOTmaskedsubs;
distribution_descriptors={'mean', 'median','std'};
output_dir3=fullfile(fMRI_threshold_dir, 'withLEXandPER_tval_1.65_threshold_posonly');
processing_algo='tavl_threshold_posonly';
tval_threshold=1.65; 
sub_num_threshold3=1:5:20;

[subs_VOTfMRI3,minmaxval3] = processingSubstruct(VOTmaskedsubs,tempmgh,tval_threshold,processing_algo);
[vertex_struct3,vertex_valuecount3] = transVertexstruct(subs_VOTfMRI3,tempmgh);

createStatsmap(vertex_struct3,tempmgh,sub_num_threshold3,distribution_descriptors, output_dir3,tval_threshold)



%% draw subject overlap histogram from 9227 vertex
inflated_brain_dir= '/export/home/tlei/public/Gari/MINI/ANALYSIS/freesurferacpc/fsaverage/surf/lh.inflated';
mgh_dir='/export/home/tlei/public/Gari/MINI/ANALYSIS/fMRI_SPM/block/analysis_block_acpc_lhITfusLatOcc/smoothed_fhmw5/S001/results';
mghFname='RWvsCB305.fhmw5.mgh' ;

cmd='freeview -v';

%% some update on pure fMRI analysis
% I can do first no threshold, then I will check if it is normaly
% distributed
activation_threshold=0;
[subs_VOTfMRI,vertex_struct,vertex_valuecount] = getVertexstruct(VOTmaskedsubs,tempmgh,activation_threshold);

[h,p,ksstat,cv]=kstest(vertex_struct(59238).VOT_block_RWvsFF);


%% Do the PCA and plot the map
myMultimodalPCA()
%% Replicate mri_glmfit 
myMRIGlmfit()

%%
% add a new function, just to store the thing into a mat so that when
% matlab crashed, I can still run?  or I can run it in python
%% Create probabilistic fMRI maps inside ROIs for LEX and PER contrasts
% It will do the analysis inside the adjLOTS, which is specified in kkvertex
myfMRICreateProbabilistic('TEST', TESTsubs, kkvertex, tempmgh, sm);
% Create the RETEST files as well for supplementary materials
myfMRICreateProbabilistic('RETEST', DAY2subs, kkvertex, tempmgh, sm);

%% Behavioral analysis
% We need to create the maps and correct them with MonteCarlos
% Then use the ROIs where the correlation is to obtain the average fMRI value
% And send it to R to do the plots
myBehavfMRIRegression()

%% DWI Analysis
% Create the probabilistic maps of the tracts
% For the 66 TEST datasets
myCreateDWIfsaverage('TEST', TESTsubs, tempmgh);
% Create the RETEST files as well for supplementary materials
myCreateDWIfsaverage('RT', DAY2subs, tempmgh);

% Create dicotomic variable (for Figure 4 analyses)
myDWIfMRIDicotomicData('TEST'  , TESTsubs, kkvertex)
myDWIfMRIDicotomicData('RETEST', DAY2subs, kkvertex)

% Compare the LEX and PER signals inside the cortical endings
myDWIfMRIcorticalEndings('TEST'  , TESTsubs)
myDWIfMRIcorticalEndings('RETEST'  , DAY2subs) % change the subnames to index or recode with tables

%% qMRI Analysis, create values for Figure 5 (in R)
% Create the differences between ROIs
[H,P,CI,STATS] = myqMRIAnalysisInROIs('TEST', TESTsubs, kkvertex, tempmgh);
% Create the RETEST files as well for supplementary materials
[H,P,CI,STATS] = myqMRIAnalysisInROIs('RETEST', DAY2subs, kkvertex, tempmgh);

%% Create the datafiles with all the pervertex spmT values for Figure 6 (in R)
TH = 0;
mySavefsavgSurfaceValues('TEST', TESTsubs, TH)
mySavefsavgSurfaceValues('RETEST', DAY2subs, TH)
%%
[m,i]=max(max([VOTmaskedsubs.('VOT_block_RWvsPS')]))
s31=VOTmaskedsubs(27).VOT_block_RWvsPS

logP = tempmgh;
logP.vol = [s31];
inFile = ['S031_RWvsPS'];

Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
MRIwrite(logP, [Resultdir fsp inFile '.mgh']);
