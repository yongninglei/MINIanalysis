% this is to investigate the distribution of vertex along the y axis
clc; clear all;
tbUse MINIanalysis
fsp = filesep;
% remember to setenv!!! 
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
labeldir      = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc/fsaverage/label';
% get the smoothing option
smmm = '5';
sm = '.fhmw5'; % '.fhmw5': Read smoothed data, '' read non smoothed data

% Read a reference mgh to be used afterwards
tempmgh = MRIread([fmriresultdir fsp 'S001' fsp 'results' fsp 'AllNoFacesvsNull305.mgh']);
tempmgh.vol = zeros(size(tempmgh.vol));


VOT= myFSread_label('fsaverage',[labeldir fsp 'lh.ITfusiLatOccNoV1V2yMin.label'], 1); 
VOT(:,1)=VOT(:,1) + 1;
LOTS=myFSread_label('fsaverage',[labeldir fsp 'lh.LOTS.label'], 1); 
LOTS(:,1) =LOTS(:,1)+ 1;

% load the data, list all the parameters that we need to use
fname='VOTmaskedsubs_CURVCTqMRIT1qMRIMTV6fMRI.mat';
VOTmaskedsubs=load(fullfile(ANALYSISdir,fname)).VOTmaskedsubs;
[vertex_struct,vertex_valuecount] = transVertexstruct(VOTmaskedsubs,tempmgh);

% first, reorder the LOTS based on y axis, then get every vertex for the violin plot
y_dec=sortrows(LOTS,3);

%% get only the positive value of the LEX and PER


perlex_struct=repmat(struct('pval',[], 'neglog10pval',[]), 163842,1);


for i = 1: length(vertex_struct)
    array_CB=vertex_struct(i).VOT_block_RWvsCB;
    array_PS=vertex_struct(i).VOT_block_RWvsPS;
    array_SD=vertex_struct(i).VOT_block_RWvsSD;
    
    
    array_CS=vertex_struct(i).VOT_block_RWvsCS;
    array_FF=vertex_struct(i).VOT_block_RWvsFF;
    array_PW=vertex_struct(i).VOT_block_RWvsPW;
    
    array_lex= [array_CS,array_FF,array_PW];
    array_per= [array_CB,array_PS,array_SD];
    
    mean_lex= mean(array_lex,2,'omitnan');
    mean_per= mean(array_per,2,'omitnan');
    
    [~,p,~,~]=ttest2(mean_per,mean_lex);
    if mean(mean_per)>0
        perlex_struct(i).pval=p;
        perlex_struct(i).neglog10pval=-log10(p);
    else
        perlex_struct(i).pval=0;
        perlex_struct(i).neglog10pval=0;
    end
end

%% change 1, remove the PW
per3lex2_struct=repmat(struct('pval',[], 'neglog10pval',[]), 163842,1);


for i = 1: length(vertex_struct)
    array_CB=vertex_struct(i).VOT_block_RWvsCB;
    array_PS=vertex_struct(i).VOT_block_RWvsPS;
    array_SD=vertex_struct(i).VOT_block_RWvsSD;
    
    
    array_CS=vertex_struct(i).VOT_block_RWvsCS;
    array_FF=vertex_struct(i).VOT_block_RWvsFF;
    
    array_lex= [array_CS,array_FF];
    array_per= [array_CB,array_PS,array_SD];
    
    mean_lex= mean(array_lex,2,'omitnan');
    mean_per= mean(array_per,2,'omitnan');
    
    [~,p,~,~]=ttest2(mean_per,mean_lex);
    if mean(mean_per)>0
        per3lex2_struct(i).pval=p;
        per3lex2_struct(i).neglog10pval=-log10(p);
    else
        per3lex2_struct(i).pval=0;
        per3lex2_struct(i).neglog10pval=0;
    end
end

% Write mgh
logP = tempmgh;
logP.vol = [per3lex2_struct.neglog10pval];
inFile = 'PER3_minus_LEX2';

Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
MRIwrite(logP, [Resultdir fsp inFile '.mgh']);
%% change 2, 3 by 3 = 9 contrasts
perxlex_struct=repmat(struct('CBvsFF',[],'CBvsCS',[],'CBvsPW',[], ...
                            'PSvsFF',[],'PSvsCS',[],'PSvsPW',[], ...
                            'SDvsFF',[],'SDvsCS',[],'SDvsPW',[]), 163842,1);


for i = 1: length(vertex_struct)
    array_CB=vertex_struct(i).VOT_block_RWvsCB;
    array_PS=vertex_struct(i).VOT_block_RWvsPS;
    array_SD=vertex_struct(i).VOT_block_RWvsSD;
    
    
    array_CS=vertex_struct(i).VOT_block_RWvsCS;
    array_FF=vertex_struct(i).VOT_block_RWvsFF;
    array_PW=vertex_struct(i).VOT_block_RWvsPW;
    
    [~,p_cbff,~,~]=ttest2(array_CB,array_FF);
    [~,p_cbcs,~,~]=ttest2(array_CB,array_CS);
    [~,p_cbpw,~,~]=ttest2(array_CB,array_PW);
    [~,p_psff,~,~]=ttest2(array_PS,array_FF);
    [~,p_pscs,~,~]=ttest2(array_PS,array_CS);
    [~,p_pspw,~,~]=ttest2(array_PS,array_PW);
    [~,p_sdff,~,~]=ttest2(array_SD,array_FF);
    [~,p_sdcs,~,~]=ttest2(array_SD,array_CS);
    [~,p_sdpw,~,~]=ttest2(array_SD,array_PW);
    
    if mean(array_CB)>0
        perxlex_struct(i).CBvsFF=-log10(p_cbff);
        perxlex_struct(i).CBvsCS=-log10(p_cbcs);
        perxlex_struct(i).CBvsPW=-log10(p_cbpw);
    else
        perxlex_struct(i).CBvsFF=0;
        perxlex_struct(i).CBvsCS=0;
        perxlex_struct(i).CBvsPW=0;
    end

    
    if mean(array_PS)>0
        perxlex_struct(i).PSvsFF=-log10(p_psff);
        perxlex_struct(i).PSvsCS=-log10(p_pscs);
        perxlex_struct(i).PSvsPW=-log10(p_pspw);
    else
        perxlex_struct(i).PSvsFF=0;
        perxlex_struct(i).PSvsCS=0;
        perxlex_struct(i).PSvsPW=0;
    end
    
    if mean(array_SD)>0
        perxlex_struct(i).SDvsFF=-log10(p_sdff);
        perxlex_struct(i).SDvsCS=-log10(p_sdcs);
        perxlex_struct(i).SDvsPW=-log10(p_sdpw);
    else
        perxlex_struct(i).SDvsFF=0;
        perxlex_struct(i).SDvsCS=0;
        perxlex_struct(i).SDvsPW=0;
    end


end

field_names= fieldnames(perxlex_struct);

for i=1:length(field_names)
    logP = tempmgh;
    logP.vol = [perxlex_struct.(field_names{i})];
    inFile = field_names{i};

    Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
    MRIwrite(logP, [Resultdir fsp inFile '.mgh']);

end
%% play with the numbers: if we do POS only and we do a unpaired t-test, what will happen?
% first, do this using the orignal thing
i=114227;
array_CB=vertex_struct(i).VOT_block_RWvsCB;
array_PS=vertex_struct(i).VOT_block_RWvsPS;
array_SD=vertex_struct(i).VOT_block_RWvsSD;


array_CS=vertex_struct(i).VOT_block_RWvsCS;
array_FF=vertex_struct(i).VOT_block_RWvsFF;
array_PW=vertex_struct(i).VOT_block_RWvsPW;

anova_array=[array_CS,array_FF,array_PW,array_CB,array_PS,array_SD];
catagory=repmat([1,1,1,2,2,2],58,1);
contrast=remmat(1:6,58,1);

[p, tbl, stats] = anovan(reshape(anova_array, [], 1), {reshape(contrast, [], 1), reshape(catagory, [], 1)}, 'model', 'interaction', 'varnames', {'treatment', 'category'});

% second, 

%% visualize the result, first do clustering and then make label
% Write mgh
logP = tempmgh;
logP.vol = [perlex_struct.neglog10pval];
inFile = 'PER_minus_LEX';

Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
MRIwrite(logP, [Resultdir fsp inFile '.mgh']);



% fshome = '/Applications/freesurfer';
fssubdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc';
% 
% setenv('FREESURFER_HOME', fshome);
setenv('SUBJECTS_DIR', fssubdir);
inFile = 'PER_minus_LEX';
Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
thp = '20';  % '13', '20'
thpcomma = '2.0';  % '1.3', '2.0'
cwpvalthresh = '0.05';  % '0.05', '0.01'
sig = 'pos';  % 'pos', 'abs', 'neg'
ver = '58sub-smooth_only_pos_per';
cd(Resultdir); if ~exist(inFile); mkdir(inFile); end

cmdsc = ['mri_surfcluster ' ...
    '--in ' inFile '.mgh ' ...
    '--csd ' fullfile('/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc/average/mult-comp-cor/fsaverage/lh/ITfusiLatOccNoV1V2yMin','fwhm05',sig,['th' thp],'mc-z.csd ') ...
    '--mask ' fullfile('/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc/average/mult-comp-cor/fsaverage/lh/ITfusiLatOccNoV1V2yMin', 'mask.mgh ') ...
    '--cwsig ' inFile fsp 'cache.th' thp '.' sig '.sig.cluster.mgh ' ...
    '--vwsig ' inFile fsp 'cache.th' thp '.' sig '.sig.voxel.mgh ' ...
    '--sum ' inFile fsp 'cache.th' thp '.' sig '.sig.cluster.summary ' ...
    '--ocn ' inFile fsp 'cache.th' thp '.' sig '.sig.ocn.mgh ' ...
    '--oannot ' inFile fsp 'cache.th' thp '.' sig '.sig.ocn.annot ' ...
    '--annot aparc ' ...
    '--csdpdf ' inFile fsp 'cache.th' thp '.' sig '.pdf.dat ' ...
    '--cwpvalthresh ' cwpvalthresh ' ' ...
    '--thmin ' thpcomma ' ' ...
    '--o ' inFile fsp 'cache.th' thp '.' sig '.sig.masked.mgh ' ...
    '--no-fixmni ' ...
    '--bonferroni 2 ' ...
    '--surf white ']
system(cmdsc)
%%
meantval_struct=repmat(struct('tval_lex',[],'tval_per',[]), 163842,1);
TH=1.65;

for i = 1: length(vertex_struct)
    array_CB=vertex_struct(i).VOT_block_RWvsCB;
    array_PS=vertex_struct(i).VOT_block_RWvsPS;
    array_SD=vertex_struct(i).VOT_block_RWvsSD;
    
    
    array_CS=vertex_struct(i).VOT_block_RWvsCS;
    array_FF=vertex_struct(i).VOT_block_RWvsFF;
    array_PW=vertex_struct(i).VOT_block_RWvsPW;
    
    array_CB(array_CB<TH)=nan;
    array_PS(array_PS<TH)=nan;
    array_SD(array_SD<TH)=nan;

    array_CS(array_CS<TH)=nan;
    array_FF(array_FF<TH)=nan;
    array_PW(array_PW<TH)=nan;

    array_lex= [array_CS,array_FF,array_PW];
    array_per= [array_CB,array_PS,array_SD];
    
    mean_lex= mean(array_lex,2,'omitnan');
    mean_per= mean(array_per,2,'omitnan');
    
    gpmean_mean_lex=mean(mean_lex,'omitnan');
    gpmean_mean_per=mean(mean_per,'omitnan');
    if isnan(gpmean_mean_lex)
        gpmean_mean_lex=0;
    end

    if isnan(gpmean_mean_per)
        gpmean_mean_per=0;
    end
    meantval_struct(i).tval_lex=gpmean_mean_lex;
    meantval_struct(i).tval_per=gpmean_mean_per;
end

% Write mgh
logP = tempmgh;
logP.vol = [meantval_struct.tval_lex];
inFile = ['T_threshold_165_mean_tval_LEX'];

Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
MRIwrite(logP, [Resultdir fsp inFile '.mgh']);

logP = tempmgh;
logP.vol = [meantval_struct.tval_per];
inFile = ['T_threshold_165_mean_tval_PER'];

Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
MRIwrite(logP, [Resultdir fsp inFile '.mgh']);

logP = MRIread([Resultdir fsp 'T_threshold_0_mean_tval_PER' '.mgh']);
logP.vol = (logP.vol)*0.3;
inFile = 'T_threshold_0_mean_tval_PER_03';

Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
MRIwrite(logP, [Resultdir fsp inFile '.mgh']);
%% increase the LEX by multiply it or reduce the per


per_reduce_lex_struct=repmat(struct('pval',[], 'neglog10pval',[]), 163842,1);


for i = 1: length(vertex_struct)
    array_CB=vertex_struct(i).VOT_block_RWvsCB;
    array_PS=vertex_struct(i).VOT_block_RWvsPS;
    array_SD=vertex_struct(i).VOT_block_RWvsSD;
    
    
    array_CS=vertex_struct(i).VOT_block_RWvsCS;
    array_FF=vertex_struct(i).VOT_block_RWvsFF;
    array_PW=vertex_struct(i).VOT_block_RWvsPW;
    
    array_lex= [array_CS,array_FF,array_PW];
    array_per= [array_CB,array_PS,array_SD];
    
    mean_lex= mean(array_lex,2,'omitnan');
    mean_per= mean(array_per,2,'omitnan');
    
    [~,p,~,~]=ttest2(mean_per*0.5,mean_lex);
    if mean(mean_per)>0
        per_reduce_lex_struct(i).pval=p;
        per_reduce_lex_struct(i).neglog10pval=-log10(p);
    else
        per_reduce_lex_struct(i).pval=0;
        per_reduce_lex_struct(i).neglog10pval=0;
    end
end

logP = tempmgh;
logP.vol = [per_reduce_lex_struct.neglog10pval];
inFile = ['0.5_reduced_PER-LEX'];

Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
MRIwrite(logP, [Resultdir fsp inFile '.mgh']);
%%
ind_contast_struct=repmat(struct('RWvsCB',[],'RWvsPS',[], 'RWvsSD',[],'RWvsCS',[],'RWvsFF',[], 'RWvsPW',[]), 163842,1);
TH=0;

for i = 1: length(vertex_struct)
    array_CB=vertex_struct(i).VOT_block_RWvsCB;
    array_PS=vertex_struct(i).VOT_block_RWvsPS;
    array_SD=vertex_struct(i).VOT_block_RWvsSD;
    
    
    array_CS=vertex_struct(i).VOT_block_RWvsCS;
    array_FF=vertex_struct(i).VOT_block_RWvsFF;
    array_PW=vertex_struct(i).VOT_block_RWvsPW;
    
    array_CB(array_CB<TH)=nan;
    array_PS(array_PS<TH)=nan;
    array_SD(array_SD<TH)=nan;

    array_CS(array_CS<TH)=nan;
    array_FF(array_FF<TH)=nan;
    array_PW(array_PW<TH)=nan;

    
    gpmean_mean_RWCB=mean(array_CB,'omitnan');
    gpmean_mean_RWPS=mean(array_PS,'omitnan');
    gpmean_mean_RWSD=mean(array_SD,'omitnan');
    
    gpmean_mean_RWCS=mean(array_CS,'omitnan');
    gpmean_mean_RWFF=mean(array_FF,'omitnan');
    gpmean_mean_RWPW=mean(array_PW,'omitnan');
    
    if isnan(gpmean_mean_RWCB)
        gpmean_mean_RWCB=0;
    end

    if isnan(gpmean_mean_RWPS)
        gpmean_mean_RWPS=0;
    end
    
    if isnan(gpmean_mean_RWSD)
        gpmean_mean_RWSD=0;
    end

    if isnan(gpmean_mean_RWCS)
        gpmean_mean_RWCS=0;
    end    
    if isnan(gpmean_mean_RWFF)
        gpmean_mean_RWFF=0;
    end

    if isnan(gpmean_mean_RWPW)
        gpmean_mean_RWPW=0;
    end    
    ind_contast_struct(i).RWvsCB=gpmean_mean_RWCB;
    ind_contast_struct(i).RWvsPS=gpmean_mean_RWPS;
    ind_contast_struct(i).RWvsSD=gpmean_mean_RWSD;
    ind_contast_struct(i).RWvsCS=gpmean_mean_RWCS;
    ind_contast_struct(i).RWvsFF=gpmean_mean_RWFF;
    ind_contast_struct(i).RWvsPW=gpmean_mean_RWPW;
end


field_names= fieldnames(ind_contast_struct);

for i=1:length(field_names)
    logP = tempmgh;
    logP.vol = [ind_contast_struct.(field_names{i})];
    inFile = ['individual_contrast_0thresh_' field_names{i}];

    Resultdir='/bcbl/home/public/Gari/MINI/ANALYSIS/fs_glmfit/per-lex';
    MRIwrite(logP, [Resultdir fsp inFile '.mgh']);

end
