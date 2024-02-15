%% Set up folders and variables
% Script based on the exploratory analysis file ~/code/MATLAB/myGlmfitPerVertex.m
% It has been reduced to contain only the code needed to reproduce the PNAS
% paper results.
% TODO1 = Use tables instead of structs... 

tbUse paper-rightMINI;
fsp = filesep;


% Folders
MINIdir = '/bcbl/home/public/Gari/MINI';
% where the data are stoted, point to a directory with a lot of subjects
DATAdir = [MINIdir fsp 'ANALYSIS' fsp "fMRI_SPM"];
% fMRI result subject dir
fmriresultdir= '/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/block/analysis_block_acpc_lhITfusLatOcc/SUBJECTS'
ANALYSISdir = [MINIdir fsp 'ANALYSIS' fsp 'tigerGLMFIT'];
fsdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc';

% remember to setenv!!! 
% those two are not needed in lmc02
%fsbin = '/Applications/freesurfer/bin';
%fshome = '/Applications/freesurfer'; 

%% what is this doing??? 
addpath(genpath('~/soft/export_fig'))
%%
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
 
%% Initialize structures and read them 
% Read all subjects in the folder and delete the non interesting fields of the
% struct
cd('/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/block/analysis_block_acpc_lhITfusLatOcc/SUBJECTS');
subs = dir('S*');
delFields  = {'date', 'bytes', 'isdir', 'datenum'};
subs = rmfield(subs, delFields);
% Read a reference mgh to be used afterwards
tempmgh = MRIread([fmriresultdir fsp 'S001' fsp 'results' fsp 'AllNoFacesvsNull305.mgh']);
tempmgh.vol = zeros(size(tempmgh.vol));

% Initialize structure with NaNs, and fill with data if it exists
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
% TESTsubs = subs(TESTind);
% DAY2subs = subs(DAY2ind);
% ALLsubs  = subs;
% subject_index = ALLind;

% right now didn't touch the QA of qMRI and DWI so I tick out the
% participants

trt= 'TESTgroup';
keep= true(size(TESTind));
keep([4,13,18,29,32,48,56,62])=false;

subject_index= TESTind(keep)
TESTsubs = subs(subject_index);

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
