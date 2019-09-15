function myCreateProbOverlay(maximas, glm, conName, fsdir, basedir, versionNum, ...
                             Tmin, tmplate305, dilateBy)
%% Funcion para escribir overlays para visualizar en individual y fsaverage space
% maximas: 1per every contrast, with all data for all subjects 
% glm: event, block...
% conName: PSvsNull, ...
% fsdir: FREESURFER_HOME
% basedir: analysis dir, until SUBJECTS folder included
% versionNum: just inc ase , version the excel has been saved
% Tmin: if T below this value, the vertex for this subject won't be counted
% tmplate: send one mgh volume in fsaverage space
% 
% TODOs:

% Dirs
fsp = filesep;
cd(fsdir)
setenv('FREESURFER_HOME', fsdir)
outputDir = [fsdir fsp 'fsaverage' fsp 'label' fsp 'MINIprobOverlay'];
if ~exist(outputDir)
    mkdir(outputDir)
end


% READ  FILES
lhwhite305 = read_surf([fsdir fsp 'fsaverage' fsp 'surf' fsp 'lh.white']);

%% Set freesurfer env variables
setenv('SUBJECTS_DIR', fsdir);
fshome = '/opt/freesurfer-5.3.0/freesurfer';
setenv('FREESURFER_HOME', fshome);

%% Do it in every ROI
% Lo de a,c,pVWFA esta en myCreateLabels.m                
% dilateLabelBy = {'4', '8', '16', '32'};  % Ademas existe el 1 con un vertex
dilateLabelBy = {'16'};
VWFAletter = {'vof', 'parc', 'votparc', 'aa', 'ca', 'cp', 'pp'};
autoprefijo = {'GM', 'GMyMin'};
for jj = 1:length(VWFAletter) 
    for kk=1:length(dilateLabelBy)
        autoprefijo =  [autoprefijo, [VWFAletter{jj} dilateLabelBy{kk}]];
    end
end


for jj = 1:length(autoprefijo) 
    prefijo = autoprefijo{jj};
    % Make the template 0
    tmplate305.vol(:) = 0;
    tmplate305.fspec  = ' ';
    tmplate305.pwd    = ' ';
    % Read the vertex values for the roi 
    ROIname = [prefijo '_vtx305' ]
    vtx_list = arrayfun(@(i) maximas{i}.(ROIname), 1:numel(maximas), 'UniformOutput', false); 
    vtxs = cell2mat(vtx_list);
    % Read the vertex values for the roi 
    TROIname = [prefijo '_T' ]
    T_list = arrayfun(@(i) maximas{i}.(TROIname), 1:numel(maximas), 'UniformOutput', false); 
    Ts = cell2mat(T_list);
    
    % Threshold values with the Tmin
    vtxsTmin = vtxs( logical([~(Ts < Tmin)] .* [~isnan(Ts)] ));
    
    
    % Count the unique ones
    ux = unique(vtxsTmin);
    if length(ux) == 1, counts = length(vtxsTmin);
    else counts = hist(vtxsTmin, ux); end
    [C] = unique(vtxsTmin', 'rows');
    toWrite = [C, counts'];
    % Create the overlay and write it
    tmplate305.vol(1, [(C+1)']) = counts/size(vtxsTmin, 2); %C+1 pq Malab es base 1
    filename = [outputDir fsp glm '_' conName '_' prefijo versionNum '.mgh'];
    MRIwrite(tmplate305, filename);
    % Create labels, dilate them, and repeat
    outputDir = [fsdir fsp 'fsaverage' fsp 'label' fsp 'MINIprobOverlay'];
    tmpDir = [outputDir fsp 'tmpLabels_' glm '_' conName '_' prefijo versionNum];
    if ~exist(tmpDir), system(['mkdir ' tmpDir]); end
    dilatedVtxs = [];
    for nl = 1:length(toWrite)
        vertex = toWrite(nl, 1);
        zenbat = toWrite(nl, 2);
        ScannerRAS = lhwhite305(vertex + 1, :); % ajustando para matlab a base 1
        subname = 'fsaverage';
        izena = ['vertex_' num2str(vertex)];
        ok = write_label(vertex, ... % Aqui no hay que ajustar, es base 0
                 ScannerRAS, ...
                 1, ...
                 [tmpDir fsp izena '.label'], ...
                 subname);
        % Y ahora lo dilato por dilateBy = 4
        system([fshome fsp 'bin' fsp 'mris_label_calc ' ...
                'dilate ' dilateBy ' ' ...
                tmpDir fsp izena ' ' ...
                fsdir fsp subname fsp 'surf' fsp 'lh.white ' ...
                tmpDir fsp izena '_dilBy' dilateBy ...
                ]);
        ROI = myFSread_label(subname, ...
                             [tmpDir fsp izena '_dilBy' dilateBy '.label'], ...
                             1);
        ROIvtx = ROI(:, 1)';
        dilatedVtxs = [dilatedVtxs, repmat(ROIvtx, 1, zenbat)];
    end
    % Count the unique ones
    dil_ux = unique(dilatedVtxs);
    if length(dil_ux) == 1, counts = length(dilatedVtxs);
    else dil_counts = hist(dilatedVtxs, dil_ux); end
    [dil_C] = unique(dilatedVtxs', 'rows');
    dil_toWrite = [dil_C, dil_counts'];
    % Make the template 0
    tmplate305.vol(:) = 0;
    % Create the overlay and write it with the dilated label data
    tmplate305.vol(1, [(dil_C+1)']) = dil_counts/size(vtxs, 2); % Matlab base 1 ver arriba
    minimo = floor(((min(dil_counts) / size(vtxs, 2))) *100);
    maximo = ceil( (max(dil_counts) / size(vtxs, 2))*100);
    
    dil_filename = [outputDir fsp glm '_' conName '_dilBy' dilateBy ...
                    '_min' sprintf('%03d', minimo) '_max' sprintf('%03d', maximo) ...
                    '_' prefijo versionNum '.mgh'];
    MRIwrite(tmplate305, dil_filename);
end
        
        
 

end

% subname = 'fsaverage'
% [status, results] = system(['freeview -viewport 3d ' ...
%         '-f ' fsdir fsp subname fsp 'surf' fsp ...
%         'lh.inflated:annot=aparc.annot' ...
%         ':overlay=' filename ' &' ...
%        ])  





