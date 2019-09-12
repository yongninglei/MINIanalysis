%% Crear ROIs reduciendo el tamano de los tractos, en superficie


% Despues de haber hecho el ROI analisis en funcional, mi objetivo en este
% caso es el:
% 1.- Encontrar las fibras del VOF y del poserior arcuate y del arcuate
% 2.- HAcer conteo de fibras, caracteristicas, y crear ROIs de los tractos
% 3.- Crear ROI-s individuales de donde estan llegando estos tractos


% Folder Names
fsp = filesep;
AnalysisDir = '/bcbl/home/public/Gari/MINI/ANALYSIS';
fs_SUBJECTS_DIR = fullfile(AnalysisDir, 'freesurferacpc');
DWIdir  = fullfile(AnalysisDir, 'DWI');
cd(DWIdir);
subs = dir('S*');
retDIR = fullfile(AnalysisDir, 'ret');
fMRIDIR = fullfile(AnalysisDir, 'fMRI_SPM');
fsbin = '/opt/freesurfer-5.3.0/freesurfer/bin';
fshome = '/opt/freesurfer-5.3.0/freesurfer'; 

MNI305to152 =     [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840];
MNI305to152sq =   [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840
                    0             0         0         1  ];
      
[vtx, labels, ct] = read_annotation([fs_SUBJECTS_DIR fsp 'fsaverage' fsp ...
                            'label' fsp 'lh.VWFAaacacppp.annot']);
lhwhite = read_surf([fs_SUBJECTS_DIR fsp 'fsaverage' fsp ...
                            'surf' fsp 'lh.white']);
ppVWFA = read_label('fsaverage', 'pp16'); 
size(ppVWFA)

codigos = unique(labels);
uniqueLabels = codigos(2:end);
myLabelCoord152 = cell(4,1);

for ii=1:length(uniqueLabels)
    myLabel{ii}     = uniqueLabels(ii);
    myLabelind      = ismember(labels, myLabel{ii});
    myLabelVtx      = vtx(myLabelind);
    myLabelCoord    = lhwhite(myLabelind,:);
    myLabelCoordvtx = lhwhite(myLabelVtx + 1,:);
    
    
    myLabelCoord152{ii} = MNI305to152 * [myLabelCoord, ones(size(myLabelCoord,1),1)]';
    
    csvwrite([fMRIDIR fsp  'fsaverage305_VWFAaacacppp_' num2str(ii) '_coords.csv'], ...
             myLabelCoord152{ii}');
%     [val, xminind] = min(myLabelCoord(:,1));
%     [val, xmaxind] = max(myLabelCoord(:,1));
%     [val, yminind] = min(myLabelCoord(:,2));
%     [val, ymaxind] = max(myLabelCoord(:,2));
%     
%     xmin{ii} = myLabelCoord(xminind,:);
%     xmax{ii} = myLabelCoord(xmaxind,:);
%     ymin{ii} = myLabelCoord(yminind,:);
%     ymax{ii} = myLabelCoord(ymaxind,:);
%     
%     xmin{ii} = MNI305to152 * [xmin{ii},1]';
%     xmax{ii} = MNI305to152 * [xmax{ii},1]';
%     ymin{ii} = MNI305to152 * [ymin{ii},1]';
%     ymax{ii} = MNI305to152 * [ymax{ii},1]';
%     
%     myLabel2plot{ii} = [xmin{ii}';  xmax{ii}'; ymin{ii}';  ymax{ii}'];
%     
end

