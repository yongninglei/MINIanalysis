%% Escribir prob overlays y labels de tractos en fsaverage space
% 
% TODOs:

% Dirs and definitions
fsp = filesep;
ANALYSISdir = '/bcbl/home/public/Gari/MINI/ANALYSIS'
fMRIdir = [ANALYSISdir fsp 'fMRI_SPM']
fsdir = [ANALYSISdir fsp 'freesurferacpc']
glmAnDir = 'analysis_event_acpc_lhITfusLatOcc';
glm = 'event';
basedir = [fMRIdir fsp glm fsp glmAnDir fsp 'SUBJECTS'];

MNI305to152 =     [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840];
MNI305to152sq =   [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840
                    0             0         0         1  ];

% READ FSAVERAGE FILES
lhwhite305 = read_surf(['/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc' ...
                      '/fsaverage/surf/lh.white']);
cd(fsdir)
setenv('FREESURFER_HOME', fsdir)
outputDir = [fsdir fsp 'fsaverage' fsp 'label'];

% READ  FILES
lhwhite305 = read_surf([fsdir fsp 'fsaverage' fsp 'surf' fsp 'lh.white']);
spmTsurf305 = MRIread([basedir fsp 'S001' fsp 'results' fsp 'RWvsPS305.mgh']); 
tmplate305 = spmTsurf305;

%% Set freesurfer env variables
setenv('SUBJECTS_DIR', fsdir);
fshome = '/opt/freesurfer-5.3.0/freesurfer';
fsbin = '/opt/freesurfer-5.3.0/freesurfer/bin';
setenv('FREESURFER_HOME', fshome);

%% Do it for every ROI
% Lo de a,c,pVWFA esta en myCreateLabels.m                
% dilateLabelBy = {'4', '8', '16', '32'};  % Ademas existe el 1 con un vertex
dilateLabelBy = {'16'};
VWFAletter = {'vof', 'parc', 'votparc'};
autoprefijo = {};
for jj = 1:length(VWFAletter) 
    for kk=1:length(dilateLabelBy)
        autoprefijo =  [autoprefijo, [VWFAletter{jj} dilateLabelBy{kk} '_305']];
    end
end



for jj = 1:length(autoprefijo) 
    prefijo = autoprefijo{jj}
    % Make the template 0
    tmplate305.vol(:) = 0;
    tmplate305.fspec  = ' ';
    tmplate305.pwd    = ' ';
    % Read the vertex values for the roi 
    
    % SEPARATE IN TEST AND RETEST TRACTS
    
    % TEST
    cd(fsdir)
    groupscurly;
    subjectsGroup{1} = [ONCE;TWICE1];
    subjectsGroup{2} = [TWICE1];
    subjectsGroup{3} = [TWICE2];
    groupName = {'TEST', 'TWICE1','RETEST'};
    for ii = 1:3
        subjects = subjectsGroup{ii};
        
        AllVtx = [];
        for ns = 1: length(subjects)
                subname = subjects{ns}
                subjectLabels = read_label(subname, prefijo);
                AllVtx = [AllVtx, subjectLabels(:,1)'];
        end

        dil_ux = unique(AllVtx);
        if length(dil_ux) == 1, counts = length(AllVtx);
        else dil_counts = hist(AllVtx, dil_ux); end
        [dil_C] = unique(AllVtx', 'rows');
        dil_toWrite = [dil_C, dil_counts'];

        % Create the overlay and write it with the dilated label data
        % tmplate305.vol(1, [(dil_C+1)']) = dil_counts/size(subjects, 1); % Matlab base 1 ver arriba
        tmplate305.vol(1, [(dil_C+1)']) = dil_counts;
        gName = groupName{ii};
        N = size(subjects, 1);
        minimo = floor(((min(dil_counts) / size(subjects, 1))) *100);
        maximo = ceil( (max(dil_counts) / size(subjects, 1))*100);

        dil_filename = [outputDir fsp ...
                        'avgTract_' prefijo ...
                        '_' gName ...
                        '_N'   sprintf('%03d', N) ...
                        '_min' sprintf('%03d', minimo) ...
                        '_max' sprintf('%03d', maximo) ...
                        '.mgh'];
        MRIwrite(tmplate305, dil_filename);
    end
end
        
        

% subname = 'fsaverage'
% [status, results] = system(['freeview -viewport 3d ' ...
%         '-f ' fsdir fsp subname fsp 'surf' fsp ...
%         'lh.inflated:annot=aparc.annot' ...
%         ':overlay=' filename ' &' ...
%        ])  





