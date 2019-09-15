%% Crear ROIs reduciendo el tamano de los tractos, en superficie


% Despues de haber hecho el ROI analisis en funcional, mi objetivo en este
% caso es el:
% 1.- Encontrar las fibras del VOF y del poserior arcuate y del arcuate
% 2.- HAcer conteo de fibras, caracteristicas, y crear ROIs de los tractos
% 3.- Crear ROI-s individuales de donde estan llegando estos tractos


% Dirs and definitions
fsp = filesep;
AnalysisDir = '/bcbl/home/public/Gari/MINI/ANALYSIS';
fs_SUBJECTS_DIR = fullfile(AnalysisDir, 'freesurferacpc');
DWIdir  = fullfile(AnalysisDir, 'DWI');
cd(DWIdir);
subs = dir('S*');
retDIR = fullfile(AnalysisDir, 'ret');
fMRIDIR = fullfile(AnalysisDir, 'fMRI_SPM', 'block', 'data');
fsbin = '/opt/freesurfer-5.3.0/freesurfer/bin';
fshome = '/opt/freesurfer-5.3.0/freesurfer'; 

MNI305to152 =     [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840];
MNI305to152sq =   [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840
                    0             0         0         1  ];
% Read 305 structures
lhwhite305 = read_surf([fs_SUBJECTS_DIR fsp 'fsaverage' fsp ...
                            'surf' fsp 'lh.white']);      
[vtx305, labels305, ct305] = read_annotation([fs_SUBJECTS_DIR fsp 'fsaverage' fsp ...
                            'label' fsp 'lh.aparc.annot']);
T1305 = MRIread(['/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc' ...
                   fsp 'fsaverage' fsp 'mri' fsp 'T1.mgz']);
% ct305.table
% ct305.struct_names
fusiform         = 7+1;
lateraloccipital = 11+1;
inferiortemporal = 9+1;

codFusi = ct305.table(fusiform, 5);
codLatOcc = ct305.table(lateraloccipital, 5);
codIT = ct305.table(inferiortemporal, 5);

myROIind = ismember(labels305, [codFusi, codLatOcc, codIT]');

% INCLUIR IT Lat Occ FF
myROIlabels305 = labels305 .* myROIind;
myROIannotName = [fs_SUBJECTS_DIR fsp 'fsaverage' fsp ...
                            'label' fsp 'lh.ITfusiLatOcc.annot'];
write_annotation(myROIannotName, vtx305, myROIlabels305, ct305);
vtxITfusiLatOcc_305 = vtx305(myROIind);                      

% INCLUIR IT Lat Occ FF NO INCLUIR V1 y V2, y poner ya lo de yMin = -40
% y hacer el resto de analisis aqui tb
yMin152  = [-30; -30; -20; 1];
yMin305  = inv(MNI305to152sq) * yMin152;
yMinSurf305 = T1305.tkrvox2ras * inv(T1305.vox2ras) * yMin305;
yMin305lhwhite = yMinSurf305(2,1);

vtxITfusiLatOcc_305 = vtx305(myROIind);
V1label = read_label('fsaverage', 'lh.V1');  % Esta sin thresholdear
V1labelVtx = V1label(:,1);
V2label = read_label('fsaverage', 'lh.V2');  % Esta sin thresholdear
V2labelVtx = V2label(:,1);

guardarV1 = ismember(vtxITfusiLatOcc_305, V1label);
guardarV2 = ismember(vtxITfusiLatOcc_305, V2label);
guardar   = guardarV1 | guardarV2;
quitar = ~guardar;
vtxITfusiLatOcc_305noV1V2 = vtxITfusiLatOcc_305(quitar);
% Ahora thresholdearlo: 
lxzy305 = lhwhite305(vtxITfusiLatOcc_305noV1V2+1,:);
lhzy_yMin_ind = find(lxzy305(:,2) < yMin305lhwhite);
lxzy305yMin = lxzy305(lhzy_yMin_ind,:);
vtxITfusiLatOcc_305noV1V2yMin = vtxITfusiLatOcc_305noV1V2(lhzy_yMin_ind);
ok = write_label(vtxITfusiLatOcc_305noV1V2yMin, ...
                 lxzy305yMin, ...
                 zeros(size(vtxITfusiLatOcc_305noV1V2yMin)), ...
                 [fs_SUBJECTS_DIR fsp 'fsaverage' fsp ...
                            'label' fsp 'lh.ITfusiLatOccNoV1V2yMin.label']);
            
% NO INCLUIR
mynotROIlabels305 = labels305 .* ~myROIind;
mynotROIannotName = [fs_SUBJECTS_DIR fsp 'fsaverage' fsp ...
                            'label' fsp 'lh.notITfusiLatOcc.annot'];
write_annotation(mynotROIannotName, vtx305, mynotROIlabels305, ct305);
vtxnotITfusiLatOcc_305 = vtx305(~myROIind);                      



%%%%%   CLUSTER PARPOOL    %%%%%%
% % myclusterLocal = parcluster('local');
% % myclusterLocal.NumWorkers
[st, re] = system('qstat -g c | grep matlab.q');
[Tok, Rem] = strtok(re);
[Tok, Rem] = strtok(Rem);
[Tok, Rem] = strtok(Rem);
[Tok, Rem] = strtok(Rem);
[available] = strtok(Rem)
parpool('ips_base', str2num(available))
%%%%% END CLUSTER PARPOOL  %%%%%%

parfor ns = 1 : length(subs)
    subname = subs(ns).name
% Tuve que cambiarle el nombre para que funcionara mejor el otro script    
     lbdir = [fs_SUBJECTS_DIR fsp subname fsp 'label'];
%     
%     [status, result] = system(['cp ' lbdir fsp 'lh.ITfusiLatOccNoV1V2yMin.label ' ...
%                                 lbdir fsp 'lhVotNoV1V2yMin16.label'])
% end
         
    
    
    
    
    setenv('FREESURFER_HOME', fshome); 
    setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR); 
%     EXTRAER ANNOTATION
cmd1 = [fsbin fsp 'mri_annotation2label ' ...
        '--subject ' subname  ' ' ...
        '--hemi lh ' ...
        '--labelbase '  lbdir fsp 'aparclh-'];
system(cmd1)
% IFG
cmd2 = [fsbin fsp 'mris_label_calc union ' ...
         lbdir fsp 'aparclh--018 ' ...
         lbdir fsp 'aparclh--019 ' ...
         lbdir fsp 'templhIFG16' ...
         ];
system(cmd2)
cmd3 = [fsbin fsp 'mris_label_calc union ' ...
         lbdir fsp 'templhIFG16 ' ...
         lbdir fsp 'aparclh--020 ' ...
         lbdir fsp 'lhIFG16' ...
         ];
system(cmd3)
%  PPC
cmd4 = [fsbin fsp 'mris_label_calc union ' ...
         lbdir fsp 'aparclh--008 ' ...
         lbdir fsp 'aparclh--031 ' ...
         lbdir fsp 'lhPPC16' ...
         ];
system(cmd4)

    
    
    
end   
    
    
    
    
    
    
for ns = 1 : length(subs)
    subname = subs(ns).name
    
%     AQUI DEBAJO DWI

    dmridir = fullfile(DWIdir, subname, 'dmri');
    fiberRois = {'L_VOF', 'L_Arcuate_Posterior', 'L_posteriorArcuate_vot'};
    for fr =1:length(fiberRois)
        fg = fiberRois{fr};
        oname      = fullfile(dmridir, [fg '_tracts.mgh']);
        oname305   = fullfile(dmridir, [fg '_tracts305.mgh']);

        % En espacio individual
        % INCLUIR TODO DENTRO DE LAT OCC IT FUSIFORM
        TalXFM = xfm_read([fs_SUBJECTS_DIR fsp subname fsp 'mri' ... 
                           fsp 'transforms' fsp 'talairach.xfm']);
        lhwhite = read_surf([fs_SUBJECTS_DIR fsp subname fsp ...
                            'surf' fsp 'lh.white']);
        T1 = MRIread([fs_SUBJECTS_DIR fsp subname fsp ...
                            'mri' fsp 'T1.mgz']);                        
        [vtx, labels, ct] = read_annotation([fs_SUBJECTS_DIR fsp subname fsp ...
                            'label' fsp 'lh.aparc.annot']);
        myROIind = ismember(labels, [codFusi, codLatOcc, codIT]');
        myROIlabels = labels .* myROIind;
        myROIannotName = [fs_SUBJECTS_DIR fsp subname fsp ...
                                  'label' fsp 'lh.ITfusiLatOcc.annot'];
        write_annotation(myROIannotName, vtx, myROIlabels, ct);
        vtxITfusiLatOcc = vtx(myROIind); 
        % Ahora escribir el tracto troceado
        readMGH = MRIread(oname);
        tmpSurf = readMGH.vol;
        tmpSurf(1, setdiff(1:size(tmpSurf,2), [(vtxITfusiLatOcc(:,1)+1)'])) = 0;
        readMGH.vol = tmpSurf;
        MRIwrite(readMGH, [dmridir fsp 'lhITfusiLatOcc_' fg '_ROI.mgh']);        

        % INCLUIR IT Lat Occ FF NO INCLUIR V1 y V2, y poner ya lo de yMin = -30
        % y hacer el resto de analisis aqui tb
        % Calculamos la yMin
        yMin  = inv(TalXFM) * yMin305; 
        yMinSurf = T1.tkrvox2ras * inv(T1.vox2ras) * yMin;
        yMinSurflhwhite = yMinSurf(2,1);
        
        vtxITfusiLatOcc = vtx(myROIind);
        V1label = read_label(subname, 'lh.V1');  % Esta sin thresholdear
        V1labelVtx = V1label(:,1);
        V2label = read_label(subname, 'lh.V2');  % Esta sin thresholdear
        V2labelVtx = V2label(:,1);

        guardarV1 = ismember(vtxITfusiLatOcc, V1label);
        guardarV2 = ismember(vtxITfusiLatOcc, V2label);
        guardar   = guardarV1 | guardarV2;
        quitar = ~guardar;
        vtxITfusiLatOcc_noV1V2 = vtxITfusiLatOcc(quitar);
        % Ahora thresholdearlo: 
        lxzy = lhwhite(vtxITfusiLatOcc_noV1V2+1,:);
        lhzy_yMin_ind = find(lxzy(:,2) < yMinSurflhwhite);
        lxzyyMin = lxzy(lhzy_yMin_ind,:);
        vtxITfusiLatOcc_noV1V2yMin = vtxITfusiLatOcc_noV1V2(lhzy_yMin_ind);
        % Escribir el label para usarlo luego
        ok = write_label(vtxITfusiLatOcc_noV1V2yMin, ...
                         lxzyyMin, ...
                         zeros(size(vtxITfusiLatOcc_noV1V2yMin)), ...
                         [fs_SUBJECTS_DIR fsp subname fsp ...
                                    'label' fsp 'lh.ITfusiLatOccNoV1V2yMin.label']);
        % Escribir el tracto resultante con la nueva restriccion sin V1 V2
        readMGH = MRIread(oname);
        tmpSurf = readMGH.vol;
        tmpSurf(1, setdiff(1:size(tmpSurf,2), ...
                [(vtxITfusiLatOcc_noV1V2(:,1)+1)'])) = 0;
        readMGH.vol = tmpSurf;
        MRIwrite(readMGH, [dmridir fsp 'lhITfusiLatOccNoV1V2_' fg '_ROI.mgh']);   
                                
                                
                                
                                
                                
                                
                                
                                
                                
        
                                
                                
        % NO INCLUIR IT FUSIFORM LAT OCC
        mynotROIlabels = labels .* ~myROIind;
        mynotROIannotName = [fs_SUBJECTS_DIR fsp subname fsp ...
                                    'label' fsp 'lh.notITfusiLatOcc.annot'];
        write_annotation(mynotROIannotName, vtx, mynotROIlabels, ct);
        vtxnotITfusiLatOcc = vtx(~myROIind);
        % Escribir el tracto resultante con la nueva restriccion sin V1 V2,
        % Esto sera posterior parietal o IFG o lo que toque
        readMGH = MRIread(oname);
        tmpSurf = readMGH.vol;
        tmpSurf(1, setdiff(1:size(tmpSurf,2), ...
                [(vtxnotITfusiLatOcc(:,1)+1)'])) = 0;
        readMGH.vol = tmpSurf;
        MRIwrite(readMGH, [dmridir fsp 'lhNotVot_' fg '_ROI.mgh']);   
                    

        
        
        
        
        
        
        
        
        
        
        % Ahora en espacio fsaverage 305
        % Escribimos en VOT
        readMGH305 = MRIread(oname305);
        tmpSurf305 = readMGH305.vol;
        tmpSurf305(1, setdiff(1:size(tmpSurf305,2), [(vtxITfusiLatOcc_305(:,1)+1)'])) = 0;
        readMGH305.vol = tmpSurf305;
        MRIwrite(readMGH305, [dmridir fsp 'lhITfusiLatOcc_' fg '_ROI_305.mgh']);
        
        % Escribimos en VOT sin V1V2
        readMGH305 = MRIread(oname305);
        tmpSurf305 = readMGH305.vol;
        tmpSurf305(1, setdiff(1:size(tmpSurf305,2), [(vtxITfusiLatOcc_305noV1V2yMin(:,1)+1)'])) = 0;
        readMGH305.vol = tmpSurf305;
        MRIwrite(readMGH305, [dmridir fsp 'lhITfusiLatOccNoV1V2_' fg '_ROI_305.mgh']);
        
        % Escribimos en notVOT
        readMGH305 = MRIread(oname305);
        tmpSurf305 = readMGH305.vol;
        tmpSurf305(1, setdiff(1:size(tmpSurf305,2), [(vtxnotITfusiLatOcc_305(:,1)+1)'])) = 0;
        readMGH305.vol = tmpSurf305;
        MRIwrite(readMGH305, [dmridir fsp 'lhNotVot_' fg '_ROI_305.mgh']);

    end    
    % Desde las labels en espacio individual, obtener la annot
%     USAGE: mris_label2annot
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
    lbdir = [fs_SUBJECTS_DIR fsp subname fsp 'label'];
    fslbdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc/fsaverage/label'
    myVWFAacpei_annotName =  [lbdir fsp 'lh.VWFAaacacppp.annot'];
    fslut = [fshome fsp 'FreeSurferColorLUT.txt'];

    VWFALABELS = {'aa16.label', 'ca16.label', 'cp16.label', 'pp16.label'};

%     for vl = 1:length(VWFALABELS)
%         cmd1 =[fsbin fsp 'mri_label2label ' ... 
%                '--srcsubject fsaverage ' ...
%                '--srclabel '  [fslbdir fsp VWFALABELS{vl}] ' ' ...
%                '--trgsubject '  subname  ' ' ...
%                '--trglabel '  [lbdir fsp VWFALABELS{vl}]   ' ' ...
%                '--hemi lh ' ...
%                '--regmethod surface' ...
%               ];
%         system(cmd1);
%     end
%     
%     cmd2 = [fsbin fsp 'mris_label2annot ' ...
%            '--s ' subname ' ' ...
%            '--h  lh ' ...
%            '--ctab  ' fslut ' ' ...
%            '--l ' [lbdir fsp 'aa16.label '] ...
%            '--l ' [lbdir fsp 'ca16.label '] ...
%            '--l ' [lbdir fsp 'cp16.label '] ...
%            '--l ' [lbdir fsp 'pp16.label '] ...
%            '--annot-path ' myVWFAacpei_annotName ...
%            ];
%        
%      system(cmd2);

% Ahora convierto los mgh de los dti a .labels.

    fiberRois = {'L_VOF', 'L_Arcuate_Posterior', 'L_posteriorArcuate_vot'};
    outputFiberRois = {'vof', 'parc', 'votparc'};
    for fr =1:length(fiberRois)
        fg = fiberRois{fr};
        outfg = outputFiberRois{fr};

        input = fullfile(dmridir, ['lhITfusiLatOcc_' fg '_ROI.mgh']);
        input305 = fullfile(dmridir, ['lhITfusiLatOcc_' fg '_ROI_305.mgh']);
        output = fullfile(lbdir, [outfg '16.label']);
        output305 = fullfile(lbdir, [outfg '16_305.label']);


        cmd3 = [fsbin fsp 'mri_cor2label ' ...
            '--i ' input ' ' ...    %  vol or surface overlay
            '--l ' output ' ' ... % labelfile % name of output file
            '--id 1 ' ...
            '--surf ' subname ' lh white ' ... % subject hemi <surf> : interpret input as surface overlay
            ];
        system(cmd3)
        
        cmd4 = [fsbin fsp 'mri_cor2label ' ...
            '--i ' input305 ' ' ...    %  vol or surface overlay
            '--l ' output305 ' ' ... % labelfile % name of output file
            '--id 1 ' ...
            '--surf fsaverage lh white ' ... % subject hemi <surf> : interpret input as surface overlay
            ];
        system(cmd4)
        
        
        % sin V1V2
        input = fullfile(dmridir, ['lhITfusiLatOccNoV1V2_' fg '_ROI.mgh']);
        input305 = fullfile(dmridir, ['lhITfusiLatOccNoV1V2_' fg '_ROI_305.mgh']);
        output = fullfile(lbdir, ['NoV1V2_' outfg '16.label']);
        output305 = fullfile(lbdir, ['NoV1V2_' outfg '16_305.label']);


        cmd5 = [fsbin fsp 'mri_cor2label ' ...
            '--i ' input ' ' ...    %  vol or surface overlay
            '--l ' output ' ' ... % labelfile % name of output file
            '--id 1 ' ...
            '--surf ' subname ' lh white ' ... % subject hemi <surf> : interpret input as surface overlay
            ];
        system(cmd5)
        
        cmd6 = [fsbin fsp 'mri_cor2label ' ...
            '--i ' input305 ' ' ...    %  vol or surface overlay
            '--l ' output305 ' ' ... % labelfile % name of output file
            '--id 1 ' ...
            '--surf fsaverage lh white ' ... % subject hemi <surf> : interpret input as surface overlay
            ];
        system(cmd6)
        
        
        %notVOT
        input = fullfile(dmridir, ['lhNotVot_' fg '_ROI.mgh']);
        input305 = fullfile(dmridir, ['lhNotVot_' fg '_ROI_305.mgh']);
        output = fullfile(lbdir, ['lhNotVot_' outfg '16.label']);
        output305 = fullfile(lbdir, ['lhNotVot_' outfg '16_305.label']);


        cmd7 = [fsbin fsp 'mri_cor2label ' ...
            '--i ' input ' ' ...    %  vol or surface overlay
            '--l ' output ' ' ... % labelfile % name of output file
            '--id 1 ' ...
            '--surf ' subname ' lh white ' ... % subject hemi <surf> : interpret input as surface overlay
            ];
        system(cmd7)
        
        cmd8 = [fsbin fsp 'mri_cor2label ' ...
            '--i ' input305 ' ' ...    %  vol or surface overlay
            '--l ' output305 ' ' ... % labelfile % name of output file
            '--id 1 ' ...
            '--surf fsaverage lh white ' ... % subject hemi <surf> : interpret input as surface overlay
            ];
        system(cmd8)
    end


end




