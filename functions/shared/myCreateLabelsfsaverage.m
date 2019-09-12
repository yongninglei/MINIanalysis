function myCreateLabelsfsaverage(TRT)
% Function for creating FS surface labels
% TRT = 'TEST'
% TRT = 'RETEST'


%% Init and folder names
clear all; close all; 
fsp = filesep;

% Folder names
fsbin = '/Applications/freesurfer/bin';
fshome = '/Applications/freesurfer';
fsdir = '~/Documents/BCBL_PROJECTS/MINI/ANALYSIS/freesurferacpc';
labeldir = fullfile(MINIPath, 'DATA', 'fslabeldir');
labeldiraparc = fullfile(labeldir,'aparcLabels');

%% READ FSAVERAGE FILES, needed for conversion
TalXFM305 = xfm_read([fsdir fsp 'fsaverage' ...
                   fsp 'mri' filesep 'transforms' ...
                   fsp 'talairach.xfm']);
T1305 = MRIread([fsdir fsp 'fsaverage' fsp 'mri' fsp 'T1.mgz']);
lhwhite305 = read_surf([fsdir fsp 'fsaverage' fsp 'surf' fsp 'lh.white']);
lhpial305 = read_surf([fsdir fsp 'fsaverage' fsp 'surf' fsp 'lh.pial']);
lhinflated305 = read_surf([fsdir fsp 'fsaverage' fsp 'surf' fsp 'lh.inflated']);
Norig305 = T1305.vox2ras;
Torig305 = T1305.tkrvox2ras;
MNI305to152 =     [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840];
MNI305to152sq =   [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840
                    0             0         0         1  ];
                
%% Create the litVWFA ROIs                
                
% Some papers report Talairach and MNI (I assume 152)
% This has been done with tal2icbm_spm (Jack Lancaster)
% tal2mni (Brett) not working, by many mm
aVWFA152 = [-45; -51; -12; 1]; % aTal = [-43, -48, -12] % AAL sub-gyral
cVWFA152 = [-45; -57; -12; 1]; % cTal = [-43, -54, -12] % AAL sub-gyral
pVWFA152 = [-45; -72; -10; 1]; % pTal = [-43, -68, -12] % AAL mid-occipital

aVWFA305 = inv(MNI305to152sq) * aVWFA152;
cVWFA305 = inv(MNI305to152sq) * cVWFA152;
pVWFA305 = inv(MNI305to152sq) * pVWFA152;


VWFA305 = {aVWFA305, cVWFA305, pVWFA305};
VWFAstring = {'aVWFA', 'cVWFA', 'pVWFA'};
dilateLabelBy = {'4'};  % , '8', '16'};

for ii = 1:length(VWFA305)
    % Create the original with only 1 vertex
    izena = VWFAstring{ii};
    ScannerRAS = VWFA305{ii}(1:3)';
    SurfaceCoord = T1305.tkrvox2ras * inv(T1305.vox2ras) * VWFA305{ii};
    Vertex = dsearchn(lhwhite305, SurfaceCoord(1:3,1)');
    subname = 'fsaverage';
    ok = write_label(Vertex - 1, ...
                     ScannerRAS, ...
                     1, ...
                     [labeldir fsp izena '1.label'], ...
                     subname);
    % And now do the dilations
    % dyld: Symbol not found: ___emutls_get_address, do it manually in the command line...     
    for jj = 1:length(dilateLabelBy)
        system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy{jj} ' ' ...
               labeldir fsp izena '1 ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               labeldir fsp izena dilateLabelBy{jj} ...
               ]);
    end
end

%% Create the fMRI GMax and mOTS and pOTS ROIs
% Read the csv created after the R analysis (file: MINI_PNAS_Analysis.Rmd)
dft = readtable(fullfile(MINIPath, 'DATA', 'fMRI', [TRT '_GMax.csv']));
dft = dft(:,{'TRT','Contrast','TYPE', 'xmean', 'ymean', 'zmean', 'xTYPE', 'yTYPE', 'zTYPE'});
dft = dft(dft.TRT==string(TRT),:);
% Create it in the same format I had before to be compatible with the code below
byContrast = dft{:,  {'xmean', 'ymean', 'zmean'}}';
byTYPE     = dft{1:2,{'xTYPE', 'yTYPE', 'zTYPE'}}'; % 1st row is PER, 2nd is LEX
coords_VOT = [byContrast,byTYPE];

pre = {'b'};  % b stands for block design. We have data for event-related as well. 
Cons = {'RWvsCB','RWvsCS','RWvsFF','RWvsPS','RWvsPW','RWvsSD','Perceptual','Lexical'};
consDesign = {};
kk = 0;
for ii = 1:length(pre)
    for jj =1:length(Cons)
        kk = kk + 1;
        conDesign{kk} = [pre{ii} '_' Cons{jj}]
    end
end
 
%% Now create one vertex based ROIs with the GMax, and dilate by 4 afterwards
pre = [TRT '_VOT'];
coords = coords_VOT;

% First convert all coordinates to MNI305
coords305 = inv(MNI305to152sq) * [coords; ones(1,length(conDesign))];

% Now create the 1 vertex label and dilate it by 4
dilateLabelBy = '4';
for ii = 1:length(conDesign)
    % Create the original with only 1 vertex
    izena = [pre '_' conDesign{ii}];
    ScannerRAS = coords305(:,ii);
    SurfaceCoord = T1305.tkrvox2ras * inv(T1305.vox2ras) * coords305(:,ii);
    Vertex = dsearchn(lhwhite305, SurfaceCoord(1:3,1)');
    subname = 'fsaverage';
    ok = write_label((Vertex - 1), ...
                     ScannerRAS(1:3)', ...
                     1, ...
                     [labeldir fsp izena '1.label'], ...
                     subname);
    disp([izena ' Vtx305: ' num2str(Vertex - 1)])          
                 
%     And now do the dilations
    setenv('FREESURFER_HOME', fshome);
    system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
               labeldir fsp izena '1 ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               labeldir fsp izena dilateLabelBy ...
               ]);
    disp (['--l ' labeldir fsp izena '1.label'])
end

































%% 

perc = [1,4,6]; % Revisar, depende de cada distinto coords_VOT
sem  = [2,3,5]; % Revisar, depende de cada distinto coords_VOT
perc = [1,3]; % Revisar, depende de cada distinto coords_VOT
sem  = [2,4]; % Revisar, depende de cada distinto coords_VOT

coords_VOT_perc = coords_VOT(:, perc);
coords_VOT_sem  = coords_VOT(:, sem);
coords_VOT_perc_avg = mean(coords_VOT_perc,2);
coords_VOT_sem_avg  = mean(coords_VOT_sem,2);

% Convert all coordinates to MNI305
coords_VOT_perc_avg_305 = inv(MNI305to152sq) * [coords_VOT_perc_avg; 1];
coords_VOT_sem_avg_305 = inv(MNI305to152sq) * [coords_VOT_sem_avg; 1];

% Obtener el vertex correspondiente
    setenv('FREESURFER_HOME', fshome);
    perc_ScannerRAS   = coords_VOT_perc_avg_305;
    sem_ScannerRAS    = coords_VOT_sem_avg_305;
    perc_SurfaceCoord = T1305.tkrvox2ras * inv(T1305.vox2ras) * coords_VOT_perc_avg_305;
    sem_SurfaceCoord  = T1305.tkrvox2ras * inv(T1305.vox2ras) * coords_VOT_sem_avg_305;
    
    perc_Vertex = dsearchn(lhwhite305, perc_SurfaceCoord(1:3,1)');
    sem_Vertex = dsearchn(lhwhite305, sem_SurfaceCoord(1:3,1)');
    
    subname = 'fsaverage';
    izena = 'lh.VOT_perc_averagesf';
    ok = write_label((perc_Vertex - 1), ...
                     perc_ScannerRAS(1:3)', ...
                     1, ...
                     [labeldir fsp izena '1f.label'], ...
                     subname);    
    dilateLabelBy = '4';
    system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
               labeldir fsp izena '1f ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               labeldir fsp izena dilateLabelBy ...
               ]);

   izena = 'lh.VOT_sem_averagesf';
   ok = write_label((sem_Vertex - 1), ...
                     sem_ScannerRAS(1:3)', ...
                     1, ...
                     [labeldir fsp izena '1f.label'], ...
                     subname);
   system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
               labeldir fsp izena '1f ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               labeldir fsp izena dilateLabelBy ...
               ]);

end
