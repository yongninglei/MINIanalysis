%% Set up folders and variables
clc; clear all;
tbUse MINIanalysis
fsp = filesep;
% remember to setenv!!! 
% setenv options
fsdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc';
setenv('SUBJECTS_DIR', fsdir);
%%
% Folders
MINIanalysis= '/bcbl/home/home_n-z/tlei/toolboxes/MINIanalysis';
MINIdir = '/bcbl/home/public/Gari/MINI';
% this is also the output dir
ANALYSISdir = [MINIdir fsp 'ANALYSIS' fsp 'individual_subject_analysis'];


%% get the subjects
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
fmri_result_dir=sprintf('/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/%s/analysis_%s_acpc_lhITfusLatOcc/SUBJECTS',design,design);
%% generate individual struct
cd(fmri_result_dir)
fsnative_subs = dir('S*');
delFields  = {'date', 'bytes', 'isdir', 'datenum'};
fsnative_subs = rmfield(fsnative_subs, delFields);

% this part need to check if it is contains all the sub groups
% # S072: qMRI missing, the pipeline did not end, usually registration problems
% # S086: qMRI missing, the pipeline did not end, usually registration problems
% # S097: el RETEST de lexical decision
% # S075:  no event related data
% the obove one need to add to the list if we want to run other analysis
subjects_to_remove={'S004', 'S013', 'S018', 'S029', 'S032', 'S048', 'S056', 'S067' };

[TESTsubs, subject_index] = getSubStruct(subjects_to_remove,fsnative_subs,TESTind);

% Load label for all the individuals "lh.S_oc-temp_lat.label"

ind_label_struct= TESTsubs;
for subIdx=1:length(ind_label_struct)
    sub=ind_label_struct(subIdx).name;
    %OTS=myFSread_label(sub, [fsdir fsp sub fsp 'label' fsp 'lh.S_oc-temp_lat.label'], 1);
    %AOS=myFSread_label(sub, [fsdir fsp sub fsp 'label' fsp 'lh.S_occipital_ant.label'], 1);
    %ITS=myFSread_label(sub, [fsdir fsp sub fsp 'label' fsp 'lh.S_temporal_inf.label'], 1);
    LOTS=myFSread_label(sub, [fsdir fsp sub fsp 'label' fsp 'lh.LOTS.label'], 1);
    
    %ind_label_struct(subIdx).OTS= OTS(:,1) + 1;
    %ind_label_struct(subIdx).AOS= AOS(:,1) + 1;
    %ind_label_struct(subIdx).ITS= ITS(:,1) + 1;
    
    ind_label_struct(subIdx).LOTS= LOTS(:,1) + 1;

    %ind_label_struct(subIdx).LVOT=[ind_label_struct(subIdx).OTS; ind_label_struct(subIdx).AOS; ind_label_struct(subIdx).ITS];
end


% initialize the struct
for ns = 1: length(TESTsubs)
    sub = TESTsubs(ns).name;
    for design = designs
        fmri_result_dir = sprintf('/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/%s/analysis_%s_acpc_lhITfusLatOcc/SUBJECTS',design{:},design{:});
        tempmgh=getTempmgh(fmri_result_dir,sub);
        NaNVector = NaN(size(tempmgh.vol'));
        %fsnative_subs(ns).CT                  = NaNVector;
        %fsnative_subs(ns).CURV                = NaNVector;
        %fsnative_subs(ns).DWI_vOF             = NaNVector;
        %fsnative_subs(ns).DWI_pARC            = NaNVector;
        %fsnative_subs(ns).DWI_vOF_label       = NaNVector;
        %fsnative_subs(ns).DWI_pARC_label      = NaNVector;
        %fsnative_subs(ns).DWI_vOF_notVotlabel = NaNVector;
        %fsnative_subs(ns).DWI_pARC_notVotlabel= NaNVector;
        %fsnative_subs(ns).qMRI_MTV            = NaNVector;
        %fsnative_subs(ns).qMRI_T1qMRI         = NaNVector;
        %fsnative_subs(ns).qMRI_MTV_WM         = NaNVector;
        %fsnative_subs(ns).qMRI_T1qMRI_WM      = NaNVector;
        for area = fMRIareas; for contrast = Contrasts
            fname = [area{:} '_' design{:} '_' contrast{:}];
            TESTsubs(ns).(fname)    = NaNVector;
        end; end; 
    end
    %fsnative_subs(ns).VOL                 = NaNVector;
end
% read the data
for ns = 1: length(TESTsubs)
    sub = TESTsubs(ns).name;
    fnames = fieldnames(TESTsubs(ns));

    for ii = 3:length(fnames)
        fname = fnames{ii};

        design = fname(5:9);
        cont = fname(11:end);

        filename = [fmri_result_dir fsp sub fsp 'results' fsp cont '.mgh'];
        if exist(filename,'file')
            temp = MRIread(filename);
            TESTsubs(ns).(fname) = temp.vol';
        else
            disp([sub ' ' fname ': file doesnt exist'])
        end

    end
end
%%
% Apply OTS mask to all filed and then save
% run them once, now we just load
% OTSmaskedsubs=maskIndSurf(TESTsubs,OTS_label_struct,'OTS',fmri_result_dir);

% Apply the LVOT mask, then we will have larger space
% LVOTmaskedsubs=maskIndSurf(TESTsubs,ind_label_struct,'LVOT',fmri_result_dir);
% save([ANALYSISdir fsp 'LVOTmaskedsubs_funconly.mat'],'LVOTmaskedsubs');

% Apply the LOTS mask, then we will have larger space
%LOTSmaskedsubs=maskIndSurf(TESTsubs,ind_label_struct,'LOTS',fmri_result_dir);

LOTSmaskedsubs=load([ANALYSISdir fsp 'LOTSmaskedsubs_funconly.mat']).LOTSmaskedsubs;
% save([ANALYSISdir fsp 'LOTSmaskedsubs_funconly.mat'],'LOTSmaskedsubs');

%% FSAVERAGE MASK
% get the fsaverage struct and use the mask two types: smoothed and
% unsmoothed, I prefer unsmoothed
% Read all subjects in the folder and delete the non interesting fields of the
% struct

% Initialize structure with NaNs, and fill with data if it exists
% Now I can read qMRI MTV, T1 functional as well as curv and CT with sm 5
% Read a reference mgh to be used afterwards
tempmgh = MRIread([fmri_result_dir fsp 'S001' fsp 'results' fsp 'AllNoFacesvsNull305.mgh']);
tempmgh.vol = zeros(size(tempmgh.vol));

fsaverage_subs = load([MINIdir fsp 'ANALYSIS' fsp 'fs_glmfit' fsp 'VOTmaskedsubs_CURVCTqMRIT1qMRIMTV6fMRI.mat']).VOTmaskedsubs;
delFields  = {'CT', 'CURV','qMRI_T1qMRI'};
fsaverage_subs = rmfield(fsaverage_subs, delFields);

OTS=myFSread_label('fsaverage', [fsdir fsp 'fsaverage' fsp 'label' fsp 'lh.S_oc-temp_lat.label'], 1);
AOS=myFSread_label('fsaverage', [fsdir fsp 'fsaverage' fsp 'label' fsp 'lh.S_occipital_ant.label'], 1);
ITS=myFSread_label('fsaverage', [fsdir fsp 'fsaverage' fsp 'label' fsp 'lh.S_temporal_inf.label'], 1);
LOTS=myFSread_label('fsaverage', [fsdir fsp 'fsaverage' fsp 'label' fsp 'lh.LOTS.label'], 1);
fsaverage_label.OTS = OTS(:,1)+1;
fsaverage_label.AOS = AOS(:,1)+1;
fsaverage_label.ITS = ITS(:,1)+1;
fsaverage_label.LOTS = LOTS(:,1)+1;
fsaverage_label.LVOT= [fsaverage_label.OTS; fsaverage_label.AOS; fsaverage_label.ITS];

LOTSmasked_fsaverage_subs=maskFsaverageSurf(fsaverage_subs,fsaverage_label.('LOTS'),tempmgh);

%% Now write the new mgh and visualize it using freesurfer
fnames = fieldnames(LOTSmasked_fsaverage_subs);
outputdir=[ANALYSISdir fsp 'LOTSmasked'];
% #################for fsaverage############
for i=1:length(LOTSmasked_fsaverage_subs)
    for ii = 3:length(fnames)
        contrast=fnames{ii}(11:end);
        batchMRIwritemgh(i,LOTSmasked_fsaverage_subs,contrast,fmri_result_dir,outputdir,'fsaverage')
    end 
end
% #################for fsnative############
for i=1:length(LOTSmaskedsubs)
    for ii = 3:length(fnames)
        contrast=fnames{ii}(11:end);
        batchMRIwritemgh(i,LOTSmaskedsubs,contrast,fmri_result_dir,outputdir,'fsnative')
    end 
end

%% generate the distribution of each contrast



field_name=fieldnames(LOTSmaskedsubs);
median_struct=struct();

for i = 1:length(LOTSmaskedsubs)
    data=LOTSmaskedsubs(i);
    median_struct(i).name=LOTSmaskedsubs(i).(field_name{1});
    for j =3:numel(field_name)
        field=field_name{j};
        contrast_data=data.(field);
        medianVal=median(contrast_data(contrast_data~=0));
        
        median_struct(i).(field)=medianVal;
    end   
end




%% new thinking, for each subject, I want to have 200/1766 top vertex and to see if they are in the mOTS and pOTS
topxvertex=200/1766;
% for each subject, get the top 10% vertex? 
field_name=fieldnames(LOTSmaskedsubs);
max15=struct();

for i = 1:length(LOTSmaskedsubs)
    data=LOTSmaskedsubs(i);
    max15(i).name=LOTSmaskedsubs(i).(field_name{1});
    for j =3:numel(field_name)
        field=field_name{j};
        contrast_data=data.(field);
        max15val=prctile(contrast_data(contrast_data~=0),85);
        
        max15(i).(field)=max15val;
    end   
end
