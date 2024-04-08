%% Set up folders and variables
clc; clear all;
tbUse MINIanalysis
fsp = filesep;

% Folders
MINIanalysis= '/bcbl/home/home_n-z/tlei/toolboxes/MINIanalysis';
MINIdir = '/bcbl/home/public/Gari/MINI';
ANALYSISdir = [MINIdir fsp 'ANALYSIS' fsp 'individual_subject_analysis'];

% remember to setenv!!! 
% setenv options
fsdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc';
setenv('SUBJECTS_DIR', fsdir);


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
ContrastConNull = {'RWvsCB','RWvsCS','RWvsFF','RWvsPS','RWvsPW','RWvsSD','RWvsNull'};
Contrasts       = ContrastConNull;

% create a struct, get individual subejct and their OTS only data

% Read all subjects in the folder and delete the non interesting fields of the
% struct
design="block";
fmriresultdir=sprintf('/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/%s/analysis_%s_acpc_lhITfusLatOcc/SUBJECTS',design,design);
cd(fmriresultdir)
subs = dir('S*');
delFields  = {'date', 'bytes', 'isdir', 'datenum'};
subs = rmfield(subs, delFields);

% Read the labels created in this project, fs305 space, analysis specific
% When reading lable files, remember to add 1 to the index 
kkvertex = myReadLabels(DATAdir, labeldir, ANALYSISdir);



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