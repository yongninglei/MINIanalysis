%% Script para crear labels al rededor de coordenadas.





clear all; close all; 
fsp = filesep;


% ANALYSISdir = '/bcbl/home/public/Gari/MINI/ANALYSIS'
% fMRIdir = [ANALYSISdir filesep 'fMRI_SPM']
% fsdir = [ANALYSISdir filesep 'freesurferacpc']
% fsbin = '/opt/freesurfer-5.3.0/freesurfer/bin';
% fshome = '/opt/freesurfer-5.3.0/freesurfer'; 


fsdir = '/Applications/freesurfer/subjects';
fsbin = '/Applications/freesurfer/bin';
fshome = '/Applications/freesurfer';

% When working in the laptop (well, always) it is a bad idea to copy the labels
% into the fsaverage/label folder, they should be stored in a project based
% folder
labeldir = fullfile(MINIPath, 'DATA', 'fslabeldir');
labeldiraparc = fullfile(labeldir,'aparcLabels');

% READ FSAVERAG FILES
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


% Ellos reportan tanto tal como mni, se supone que 152
% Esto esta hecho con tal2icbm_spm de Jack Lancaster
% Si se usa tal2mni de Brett no funciona, por muchos mm
aVWFA152 = [-45; -51; -12; 1]; % aTal = [-43, -48, -12] % AAL sub-gyral
cVWFA152 = [-45; -57; -12; 1]; % cTal = [-43, -54, -12] % AAL sub-gyral
pVWFA152 = [-45; -72; -10; 1]; % pTal = [-43, -68, -12] % AAL mid-occipital

aVWFA305 = inv(MNI305to152sq) * aVWFA152;
cVWFA305 = inv(MNI305to152sq) * cVWFA152;
pVWFA305 = inv(MNI305to152sq) * pVWFA152;


VWFA305 = {aVWFA305, cVWFA305, pVWFA305};
VWFAstring = {'aVWFA', 'cVWFA', 'pVWFA'};
% VWFA305Vertex = {aVWFA305_Vertex, cVWFA305_Vertex, pVWFA305_Vertex};
dilateLabelBy = {'4', '8', '16', '32'};

for ii = 1:length(VWFA305)
    % Create the original with only 1 vertex
    izena = VWFAstring{ii};
    ScannerRAS = VWFA305{ii}(1:3)';
    SurfaceCoord = T1305.tkrvox2ras * inv(T1305.vox2ras) *VWFA305{ii};
    Vertex = dsearchn(lhwhite305, SurfaceCoord(1:3,1)');
    subname = 'fsaverage';
    ok = write_label(Vertex - 1, ...
                     ScannerRAS, ...
                     1, ...
                     [labeldir fsp izena '1.label'], ...
                     subname);
    % And now do the dilations
    % mris_label_calc dilate 8 pVWFA1  ../surf/lh.white pVWFA8
    for jj = 1:length(dilateLabelBy)
        system(['mris_label_calc dilate ' dilateLabelBy{jj} ' ' ...
               labeldir fsp izena '1 ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               labeldir fsp izena dilateLabelBy{jj} ...
               ]);
    end
end

% AHORA HACERLO CON LOS RESULTADOS DE R, mean contrasts and clusters
% Estas son las coordenadas resultado de R



% Estos son valores en horiz y con las medias the cluster incluidas, las guardo
% pero no las uso.
% BLOCK = [ 
%     -38.75823	-43.64993	-40.41244	-39.79494	-41.64231	-40.27259	-39.60242	-41.86057
% 	-71.99183	-58.80453	-58.18042	-72.52421	-62.13687	-70.25624	-71.58759	-59.57539
% 	-7.881664	-10.981442	-9.797886	-8.694976	-5.661627	-6.709939	-7.757352	-8.951527
% ];
% EVENT = [
%     -40.23103	-43.67397	-44.54275	-39.59791	-41.30531	-41.07424	-40.31793	-43.24428
% 	  -72.909	-51.38827	-48.77277	-74.76151	-57.85768	-71.05127	-72.87011	-52.46902
% 	  -7.206417	-12.340435	-12.562641	-7.445282	-8.29031	-8.045638	-7.582563	-11.122369
% ];
% coords = [BLOCK, EVENT];

pre = {'b','e'};
% Cons = {'RWvsCB','RWvsCS','RWvsFF','RWvsPS','RWvsPW','RWvsSD','Perceptual','Semantic'};
Cons = {'RWvsCB','RWvsCS','RWvsFF','RWvsPS','RWvsPW','RWvsSD'};
% Cons = {'Perceptual','Semantic'};
consDesign = {};
kk = 0;
for ii = 1:2
    for jj =1:length(Cons)
        kk = kk + 1;
        conDesign{kk} = [pre{ii} '_' Cons{jj}]
    end
end



% VOT
coords_VOT= [
  -38.7582  -71.9918   -7.8817
  -43.6499  -58.8045  -10.9814
  -40.4124  -58.1804   -9.7979
  -39.7949  -72.5242   -8.6950
  -41.6423  -62.1369   -5.6616
  -40.2726  -70.2562   -6.7099
  -40.2310  -72.9090   -7.2064
  -43.6740  -51.3883  -12.3404
  -44.5427  -48.7728  -12.5626
  -39.5979  -74.7615   -7.4453
  -41.3053  -57.8577   -8.2903
  -41.0742  -71.0513   -8.0456]'


coords_VOT= [
  -38.7582  -71.9918   -7.8817
  -43.6740  -51.3883  -12.3404
  -40.4124  -58.1804   -9.7979
  -39.7949  -72.5242   -8.6950
  -41.6423  -62.1369   -5.6616
  -40.2726  -70.2562   -6.7099
  -40.2310  -72.9090   -7.2064
  -43.6499  -58.8045  -10.9814
  -44.5427  -48.7728  -12.5626
  -39.5979  -74.7615   -7.4453
  -41.3053  -57.8577   -8.2903
  -41.0742  -71.0513   -8.0456]'


% coords_VOT= [
%     -39.60242	-71.58759 -7.757352
%     -41.86057   -59.57539 -8.951527
%     -40.31793   -72.87011 -7.582563
%     -43.24428   -52.46902 -11.122369]'




% PPC
coords_pPC = [
-45.59590 -51.17908 27.12163
-48.24867 -49.31895 27.95466
-48.30838 -46.54324 29.33624
-46.23607 -50.60383 28.15008
-44.28183 -56.85215 31.20541
-46.12907 -48.53464 28.36631
-48.10999 -45.48705 28.19076
-46.76721 -50.23414 30.96927
-49.74184 -42.65811 28.66107
-46.74192 -47.13124 28.02282
-45.10813 -54.84170 30.12941
-44.33728 -51.48814 28.83232]';

coords_pPC_pAF = [
-44.86215 -53.85233 28.63513
-46.52165 -51.71620 29.53077
-45.25310 -52.94319 30.42537
-44.74401 -53.87247 28.31015
-44.43383 -57.17376 30.17771
-43.26104 -55.04227 28.60914
-44.86761 -53.14989 27.65959
-45.58561 -52.44284 30.57345
-45.77661 -50.85401 29.64870
-43.90003 -53.13576 29.60716
-43.03109 -56.59386 30.14644
-44.16726 -53.52453 30.07381 
]';



coords_pPC_vOF = [
-22.37814  -76.84679 28.16331
-27.43828  -74.20746 29.54643
-23.11383  -75.75960 29.63735
-23.81059  -74.01663 28.01886
-26.00578  -76.65812 30.28211
-22.84474  -75.37609 29.42279
-24.74699  -73.89958 28.63548
-27.43723  -74.69108 30.06153
-25.79525  -71.61685 29.61143
-25.79814  -73.28056 28.28717
-26.43870  -76.95410 27.75228
-24.68300  -72.85234 28.51542
]';

% IFG
coords_IFG = [
-45.09436 24.77897  8.183308 
-45.70357 24.44496  8.510858 
-45.75259 24.94429 10.230478 
-45.51964 24.24195  9.581508 
-43.60386 27.63196  5.943976 
-45.60038 24.90873 10.578524 
-45.27744 23.79199 10.098160 
-46.05259 22.43957 10.283963 
-45.58692 24.18601  8.460057 
-46.14542 23.67219  9.702549 
-43.57243 26.45423  6.452946 
-45.32129 23.69132  9.139822]';


%% Crear annotations para visualizacion
% AQUI SEPARAR vOT, pPC, IFG
pre = 'VOT';
coords = coords_VOT;

% First convert all coordinates to MNI305
coords305 = inv(MNI305to152sq) * [coords; ones(1,length(conDesign))];

% Now create the 1 vertex label and dilate it by 2
dilateLabelBy = '1';
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
                     [labeldir fsp izena '1f.label'], ...
                     subname);
    disp([izena ' Vtx305: ' num2str(Vertex - 1)])          
                 
%     And now do the dilations
%     mris_label_calc dilate 8 pVWFA1  ../surf/lh.white pVWFA8
    setenv('FREESURFER_HOME', fshome);
%     system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
%                labeldir fsp izena '1 ' ...
%                fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
%                labeldir fsp izena dilateLabelBy ...
%                ]);
    disp (['--l ' labeldir fsp izena '1.label'])
end

% And now create the annotations out of the labels
% USAGE: mris_label2annot 
% 
%    --s subject : FreeSurfer subject 
%    --h hemi : hemisphere (lh or rh)
%    --ctab ctabfile : colortable (like FreeSurferColorLUT.txt)
%    --l label1 <--l label 2 ...> : label file(s)
%    --a annotname : output annotation file (hemi.annotname.annot)
%    --annot-path annotpath : full name/path of annotation file
%    --ldir labeldir : when not using --l
%    --ldir-default : use subject/labels as labeldir
%    --no-unknown : do not map unhit labels to index 0
%    --thresh thresh : threshold label by stats field
%    --maxstatwinner : keep label with highest 'stat' value
% 
%    --debug     turn on debugging
%    --noverbose turn off overlap and stat override messages
%    --checkopts don't run anything, just check options and exit
%    --help      print out information on how to use this program
%    --version   print out version and exit
% pre = 'pPC_pAF';
dilateLabelBy = '1';
setenv('FREESURFER_HOME', fshome);
    system([fsbin fsp 'mris_label2annot ' ...
            '--s fsaverage ' ...
            '--h lh ' ...
            '--ctab /Applications/freesurfer/FreeSurferColorLUT.txt ' ...
            '--a ' pre '_ContrastClusterAvgs' dilateLabelBy ' ' ...
'--l ' labeldir fsp pre '_b_RWvsCB' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_b_RWvsCS' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_b_RWvsFF' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_b_RWvsPS' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_b_RWvsPW' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_b_RWvsSD' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_e_RWvsCB' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_e_RWvsCS' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_e_RWvsFF' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_e_RWvsPS' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_e_RWvsPW' dilateLabelBy '.label ' ...
'--l ' labeldir fsp pre '_e_RWvsSD' dilateLabelBy '.label']);
    
system(['cp /Users/glerma/Documents/BCBL_PROJECTS/MINI/ANALYSIS/freesurferacpc/fsaverage/label/lh.VOT_ContrastClusterAvgs1.annot ' ...
        labeldir])
%% Ahora vamos a crear 6 ROIS, dos en cada VOT, PPC, IFG, de 1-3-8 dilated

if (1) % VOT
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

if (0) % PPC
perc = [1,4,6,7,10,12];
sem = [2,3,5,8,9,11];
coords_PPC_perc = coords_pPC_vOF(:,perc);
coords_PPC_sem  = coords_pPC_pAF(:,sem);
coords_PPC_perc_avg = mean(coords_PPC_perc,2);
coords_PPC_sem_avg  = mean(coords_PPC_sem,2);

% Convert all coordinates to MNI305
coords_PPC_perc_avg_305 = inv(MNI305to152sq) * [coords_PPC_perc_avg; 1];
coords_PPC_sem_avg_305 = inv(MNI305to152sq) * [coords_PPC_sem_avg; 1];

% Obtener el vertex correspondiente
    setenv('FREESURFER_HOME', fshome);
    perc_ScannerRAS   = coords_PPC_perc_avg_305;
    sem_ScannerRAS    = coords_PPC_sem_avg_305;
    perc_SurfaceCoord = T1305.tkrvox2ras * inv(T1305.vox2ras) * coords_PPC_perc_avg_305;
    sem_SurfaceCoord  = T1305.tkrvox2ras * inv(T1305.vox2ras) * coords_PPC_sem_avg_305;
    
    perc_Vertex = dsearchn(lhwhite305, perc_SurfaceCoord(1:3,1)');
    sem_Vertex = dsearchn(lhwhite305, sem_SurfaceCoord(1:3,1)');
    
    subname = 'fsaverage';
    izena = 'lh.PPC_perc_averages';
    ok = write_label((perc_Vertex - 1), ...
                     perc_ScannerRAS(1:3)', ...
                     1, ...
                     [fsdir fsp subname fsp 'label' fsp izena '1.label'], ...
                     subname);    
    dilateLabelBy = '3';
    system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
               fsdir fsp subname fsp 'label' fsp izena '1 ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               fsdir fsp subname fsp 'label' fsp izena dilateLabelBy ...
               ]);

   izena = 'lh.PPC_sem_averages';
   ok = write_label((sem_Vertex - 1), ...
                     sem_ScannerRAS(1:3)', ...
                     1, ...
                     [fsdir fsp subname fsp 'label' fsp izena '1.label'], ...
                     subname);
   system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
               fsdir fsp subname fsp 'label' fsp izena '1 ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               fsdir fsp subname fsp 'label' fsp izena dilateLabelBy ...
               ]);
end

if (0) % IFG
perc = [1,4,6,7,10,12];
sem = [5,11];
coords_IFG_perc = coords_IFG(:,perc);
coords_IFG_sem  = coords_IFG(:,sem);
coords_IFG_perc_avg = mean(coords_IFG_perc,2);
coords_IFG_sem_avg  = mean(coords_IFG_sem,2);

% Convert all coordinates to MNI305
coords_IFG_perc_avg_305 = inv(MNI305to152sq) * [coords_IFG_perc_avg; 1];
coords_IFG_sem_avg_305  = inv(MNI305to152sq) * [coords_IFG_sem_avg; 1];

% Obtener el vertex correspondiente
    setenv('FREESURFER_HOME', fshome);
    perc_ScannerRAS   = coords_IFG_perc_avg_305;
    sem_ScannerRAS    = coords_IFG_sem_avg_305;
    perc_SurfaceCoord = T1305.tkrvox2ras * inv(T1305.vox2ras) * coords_IFG_perc_avg_305;
    sem_SurfaceCoord  = T1305.tkrvox2ras * inv(T1305.vox2ras) * coords_IFG_sem_avg_305;
    
    perc_Vertex = dsearchn(lhwhite305, perc_SurfaceCoord(1:3,1)');
    sem_Vertex = dsearchn(lhwhite305, sem_SurfaceCoord(1:3,1)');
    
    subname = 'fsaverage';
    izena = 'lh.IFG_perc_averages';
    ok = write_label((perc_Vertex - 1), ...
                     perc_ScannerRAS(1:3)', ...
                     1, ...
                     [fsdir fsp subname fsp 'label' fsp izena '1.label'], ...
                     subname);    
    dilateLabelBy = '3';
    system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
               fsdir fsp subname fsp 'label' fsp izena '1 ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               fsdir fsp subname fsp 'label' fsp izena dilateLabelBy ...
               ]);

   izena = 'lh.IFG_sem_averages';
   ok = write_label((sem_Vertex - 1), ...
                     sem_ScannerRAS(1:3)', ...
                     1, ...
                     [fsdir fsp subname fsp 'label' fsp izena '1.label'], ...
                     subname);
   system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
               fsdir fsp subname fsp 'label' fsp izena '1 ' ...
               fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
               fsdir fsp subname fsp 'label' fsp izena dilateLabelBy ...
               ]);           
end        
 

%% LO MISMO PERO EN LOS VtxMax de la prediccion
% % Apuntar los vertex maximos
% IFG
maxVtx.WHzRT---IFG_block_RWvsCB = 20524 + 1;
maxVtx.WHzRT---IFG_block_RWvsPW = 125095 + 1;
% PPC
maxVtx.WHzRT---PPC_block_RWvsCB = 107027 + 1;
maxVtx.WHzRT---PPC_block_RWvsPW = 78310 + 1;
% VOT
maxVtx.WHzRT---VOT_block_RWvsCB = 13255 + 1;
% PW no sobrevive ningun cluster, asi que veo cual es el de sig max
% A = MRIread('/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/myGLMFIT/behav_fmri/WHzRT---VOT_block_RWvsPW-linear.fhmw5.mgh');
% [Y,I] = max(A.vol);
% I = 114365
maxVtx.WHzRT---VOT_block_RWvsPW = 114364 + 1;

% Escribir el label de 1 y de 3 en el vertex
setenv('FREESURFER_HOME', fshome);
subname = 'fsaverage';
dilateLabelBy = '3';
mvnames = fieldnames(maxVtx);
for mv = 1:length(mvnames)
    if mvnames{mv}(end-1:end) == 'PW'
        type = '_sem_';
    else
        type = '_perc_';
    end
    izena      = ['lh.ldmaxVtx' type mvnames{mv}];
    Vertex     = maxVtx.(mvnames{mv});
    ScannerRAS = lhwhite305(Vertex,:);  % creo que me he saltado un paso pero da igual, se fija en el vertex de todas maneras
    ok         = write_label((Vertex - 1), ...
                              ScannerRAS(1:3), ...
                              1, ...
                              [fsdir fsp subname fsp 'label' fsp izena '1.label'], ...
                              subname);    

    system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
           fsdir fsp subname fsp 'label' fsp izena '1 ' ...
           fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
           fsdir fsp subname fsp 'label' fsp izena dilateLabelBy ...
           ]);
end
        
% Ademas escribir la annotation para verlos todos juntos           
setenv('FREESURFER_HOME', fshome);
pre = 'lh.ldmaxVtx'
dilateLabelBy = '1';
system([fsbin fsp 'mris_label2annot ' ...
            '--s fsaverage ' ...
            '--h lh ' ...
            '--ctab /Applications/freesurfer/FreeSurferColorLUT.txt ' ...
            '--a ' pre '_ROIs_' dilateLabelBy ' ' ...
            '--l ' labeldir fsp pre '_perc_WHzRT_IFG_block_RWvsCB' dilateLabelBy '.label ' ...
            '--l ' labeldir fsp pre '_sem_WHzRT_IFG_block_RWvsPW' dilateLabelBy '.label ' ...
            '--l ' labeldir fsp pre '_perc_WHzRT_PPC_block_RWvsCB' dilateLabelBy '.label ' ...
            '--l ' labeldir fsp pre '_sem_WHzRT_PPC_block_RWvsPW' dilateLabelBy '.label ' ...
            '--l ' labeldir fsp pre '_perc_WHzRT_VOT_block_RWvsCB' dilateLabelBy '.label ' ...
            '--l ' labeldir fsp pre '_sem_WHzRT_VOT_block_RWvsPW' dilateLabelBy '.label ']);

               
%% Crear label desde listado de vertex

if (1) % IFG


% Obtener el vertex correspondiente
    setenv('FREESURFER_HOME', fshome);
%     listaVtx = [130302,64677,13879,44601,38938,112166,126277,112154,...
%                     38930,157743,87074,38916,27776,126234,157647,5184, ...
%                     142221,35403,49144,58751,78355,95161,122029,146850,...
%                     146855,146861,101067,8333,146956,122078,122074,131822];
   % Este es un nuevo path creado despues de revisar el anterior
%    listaVtx = [13879,38938,112166,126277,112154, ... 
%                     38930,  9562, 157732, 31340,135528, ...
%                     51990,  6845, 43567 ,15245 ,33697 , ...
%                     21305,157659,112092,142233,94329 , ...
%                     49161, 95161, 122029,146850,146861,  ...
%                     8333,146956,122078,122074,131822];
   listaVtx = [13879,38938,112166,126277,112154, ... 
                    38930,  9562, 157732, 31340,135528, ...
                    51990,  6845, 43567 ,15245 ,33697 , ...
                    21305,157659,112092,142233,94329]
    ScannerRAS = lhwhite305(listaVtx,:);  % creo que me he saltado un paso pero da igual, se fija en el vertex de todas maneras
    
    subname = 'fsaverage';
    izena = 'lh_listadoPPC_20vtx_V03';
    ok = write_label(listaVtx', ... % No hay qye quitar, son de freeview
                     ScannerRAS, ...
                     ones(size(listaVtx')), ...
                     [fsdir fsp subname fsp 'label' fsp izena '1.label'], ...
                     subname);    
%     dilateLabelBy = '3';
%     system([fsbin fsp 'mris_label_calc dilate ' dilateLabelBy ' ' ...
%                fsdir fsp subname fsp 'label' fsp izena '1 ' ...
%                fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
%                fsdir fsp subname fsp 'label' fsp izena dilateLabelBy ...
%                ]);

   
end        