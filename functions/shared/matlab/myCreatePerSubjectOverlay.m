function myCreateProbOverlay(glm, fsdir, basedir, versionNum, ...
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
Tmin = '2.34';
dilateBy = '4';
GLMs =       {'block', ...
              'event'};
glmAnDirs =  {'analysis_block_acpc_lhITfusLatOcc', ...
              'analysis_event_acpc_lhITfusLatOcc'}       
Tmin = '2.34';  
versionNum = 'V10'

fsdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc';
[vertices305, labels305, colortable305]  =  read_annotation([fsdir fsp ...
     'fsaverage' fsp 'label' fsp 'lh.aparc.annot']);
csvdir = [fsdir fsp 'fsaverage' fsp 'label' fsp 'MINIperSUBJECT' fsp ...
          'Tmin' Tmin '_' versionNum  fsp 'SUBJECTS'];
annotDir = [fsdir fsp 'fsaverage' fsp 'label' fsp ... 
            'MINIperSUBJECT' fsp 'Tmin' Tmin '_' versionNum fsp 'annots'];
cd(csvdir)
fshome = '/opt/freesurfer-5.3.0/freesurfer';
setenv('FREESURFER_HOME', fshome)
setenv('SUBJECTS_DIR', fsdir)


% Create ROI interest names
dilateLabelBy = {'i16', 'e16'};
VWFAletter = {'a', 'c', 'p'};
autoprefijo = {'GM', 'GMyMin'};
for jj = 1:length(VWFAletter) 
    for kk=1:length(dilateLabelBy)
        autoprefijo =  [autoprefijo, [VWFAletter{jj} dilateLabelBy{kk}]];
    end
end
% Add the _vtx and _vtx305 at the end
endings    = {'_vtx'};
endings305 = {'_vtx305'};
ROInames = {};
ROInames305 = {};
for kk = 1:length(endings)
    for mm = 1:length(autoprefijo)
        ROInames = [ROInames [autoprefijo{mm} endings{kk}]];
    end
end

for kk = 1:length(endings305)
    for mm = 1:length(autoprefijo)
        ROInames305 = [ROInames305 [autoprefijo{mm} endings305{kk}]];
    end
end


for noglm=1:length(GLMs)
    %%%%%%%%%%%%%%%%%
    glm      = GLMs{noglm}
    % glmAnDir = glmAnDirs{noglm};
    % basedir = [fMRIdir fsp glm fsp glmAnDir fsp 'SUBJECTS'];
    %%%%%%%%%%%%%%%%%
    
    

    % Read all data
    datos = dir([glm  '_*.csv']);
    for ns = 1:length(datos)
        csvName = datos(ns).name;
        subname = csvName(7:10);
        
        % Read each subjects data
        A = readtable(csvName);
        data = table2dataset(A);
    
         
        % Read aparc to substitute it with our data
        [vertices, labels, colortable]  =  read_annotation([fsdir fsp ...
                                subname fsp 'label' fsp 'lh.aparc.annot']);
                            
        % We will do: 
        % Create annotations per each subject
        % We will have, per each subject, the vtx and the vtx305 versions
        % each annotation will have N amount of contrasts per one of the M rois of
        % interest.
        % each ROI will start with a 10,20,30,40...
        % each  conName will be : 01, 02, 03
        % Then 1001 will be ROI = GM  and conName = RWvsNull, for example
        % 1000 will be the ROI itself in order to visualize it. 
        
        
        
        
        % UPDATE VERTICES
        % We don't need to do anything here
        
        
        % UPDATE COLORTABLE
  
        % Anado  y MOreThan1hit al final    
     
        structNames = {};
        for rn = 1 : length(ROInames)
              for cn = 1:length(data.Contrast)
                structNames = [structNames [ROInames{rn} '_'  data.Contrast{cn}]];
              end
        end    
        structNames305 = {};
        for rn = 1 : length(ROInames305)
              for cn = 1:length(data.Contrast)
                structNames305 = [structNames305 [ROInames305{rn} '_'  data.Contrast{cn}]];
              end
        end          
        colortable.struct_names = [structNames'; 'MoreThan1Hit'];
        colortable305.struct_names = [structNames305'; 'MoreThan1Hit'];
        
        colorsPerROI =  colortable.table(1:length(data.Contrast), 1:3);
        colorsAll    =  repmat(colorsPerROI, [length(ROInames), 1]);
        % Parece ser q el ultimo es opacity, por si acaso lo mantengo en 0
        jitter       =  rectpulse(0:length(ROInames)-1, length(data.Contrast))';
        jitter3      =  repmat(jitter, [1,3,]);
        colorsAllJitter = colorsAll + jitter3;
        tabla4       =  rectpulse(zeros(1,length(ROInames)), length(data.Contrast))';
        tabla5       =  colorsAllJitter(:,1) + (colorsAllJitter(:,2) .* 2^8) + ...
                        (colorsAllJitter(:,3) .* 2*16) + (tabla4 .* 2^24);
        % add last color for qhen more than one contrast are in the same
        % vertex
        colortable.table = [colorsAllJitter, tabla4, tabla5; 71,255,71,0,99999];
        
        colortable305.table = colortable.table;
        
        colortable.numEntries    = length(structNames)    + 1;
        colortable305.numEntries = length(structNames305) + 1;
       
        colortable.orig_tab      = [annotDir fsp glm  '_' subname '.annot.ctab'];
        colortable305.orig_tab   = [annotDir fsp glm  '_' subname '_305.annot.ctab'];
        
        
        % UPDATE LABELS
        % indiv
        labels    = zeros(size(labels));
        
        vtxAll = [];
        for rn = 1 : length(ROInames)
            vtxAll = [vtxAll; data.(ROInames{rn})];
        end
        
        vtxAll1 = [vtxAll + 1; NaN];
        
        
        matcheo = [vtxAll1, colortable.table(:,5)];
        mathNonan = matcheo(~isnan(vtxAll1), :);
        
        ux = unique(mathNonan(:,1));
        if length(ux) == 1, counts = length(mathNonan(:,1));
        else counts = hist(mathNonan(:,1), ux); end
        [C] = unique(mathNonan(:,1)', 'rows');
        
        match99999 = mathNonan;
        match99999(counts' > 1 ,2) = 99999; 

        
        labels([match99999(:,1)]) = match99999(:,2);
        
        
        
        % 305
        labels305 = zeros(size(labels305));
        vtxAll = [];
        for rn = 1 : length(ROInames305)
            vtxAll = [vtxAll; data.(ROInames305{rn})];
        end
        
        vtxAll1 = [vtxAll + 1; NaN];
        
        
        matcheo = [vtxAll1, colortable305.table(:,5)];
        mathNonan = matcheo(~isnan(vtxAll1), :);
        
        ux = unique(mathNonan(:,1));
        if length(ux) == 1, counts = length(mathNonan(:,1));
        else counts = hist(mathNonan(:,1), ux); end
        [C] = unique(mathNonan(:,1)', 'rows');
      
        
        match99999 = mathNonan;
        match99999(counts' > 1 ,2) = 99999; 

        
        labels305([match99999(:,1)]) = match99999(:,2);
        
        
        
        
        % WRITE ANNOTATION
        
        filenameAnnot = [annotDir fsp glm  '_' subname '.annot'];
        filenameAnnot305 = [annotDir fsp glm  '_' subname '_305.annot'];
        
                
        write_annotation(filenameAnnot, vertices, labels, colortable);
        write_annotation(filenameAnnot305, vertices305, labels305, colortable305);
            
            
            
   end
       
end
end

% testing
% [ver, lab, ct] = read_annotation('block_S001_305.annot');





