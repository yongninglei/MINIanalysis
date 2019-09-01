%% Script para obtener local y global maximas
% TODOs:
% OK 1.- Eliminar usar SPM para leer las GM, se puede leer del contraste,
%        asegurar que el resultado es el mismo.
% OK 2.- En vez de usar Jobard, usar solo local, pero basado en el paper de
%        Dehahene traer el anterior posterior y central/classical VWFA
% OK 3.- Dar la vuelta a :
%     3.1.- FacesvsWords, multiplicar por -1
%     3.2.- WpsPW, multiplicar por -1, queremos PW encima
%     3.3.- FAces vs PW
% 4.- Voy a tener que predecir con la tractografia la ubicacion del VWFA y
% ver con que contraste soy capaz de predecir mejor uno u otro
% 5.- Ademas hay que hacer el modelo de Brian en LiFE, cual adivina mejor
% la ubicacion, un sujeto nuevo o el modelo?
% 6.- Preguntar a Kepa aquello que me pregunto sobre lo del IFG
% OK 7.- Escribir script myCreateAnnotation.m para poder escribir en
% individual space o en fsaverage los datos y se vea como un heat map. ---
% al final he hecho overap, si no tengo que hacer un LUT, no seria un mapa
% probabilistico, con Overlay si lo puedo hacer





clear all; close all; 


% Dirs
fsp = filesep;
ANALYSISdir = '/bcbl/home/public/Gari/MINI/ANALYSIS'
fMRIdir = [ANALYSISdir fsp 'fMRI_SPM']
fsdir = [ANALYSISdir fsp 'freesurferacpc']




% READ FSAVERAGE FILES
TalXFM305 = xfm_read(['/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc' ...
                   fsp 'fsaverage' fsp 'mri' fsp 'transforms' fsp 'talairach.xfm']);
T1305 = MRIread(['/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc' ...
                   fsp 'fsaverage' fsp 'mri' fsp 'T1.mgz']);
lhwhite305 = read_surf(['/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc' ...
                      '/fsaverage/surf/lh.white']);
lhpial305 = read_surf(['/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc' ...
                      '/fsaverage/surf/lh.pial']);
lhinflated305 = read_surf(['/bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc' ...
                      '/fsaverage/surf/lh.inflated']);
Norig305 = T1305.vox2ras;
Torig305 = T1305.tkrvox2ras;
MNI305to152 =     [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840];
MNI305to152sq =   [  0.9975   -0.0073    0.0176   -0.0429
                     0.0146    1.0009   -0.0024    1.5496
                    -0.0130   -0.0093    0.9971    1.1840
                    0             0         0         1  ];
                


% Ademas buscar la y minima para que no se vayan tan anteriores
yMin152  = [-30; -40; -20; 1];
yMin305  = inv(MNI305to152sq) * yMin152;
yMinSurf305 = T1305.tkrvox2ras * inv(T1305.vox2ras) * yMin305;


%%%%%   CLUSTER PARPOOL    %%%%%%
% % myclusterLocal = parcluster('local');
% % myclusterLocal.NumWorkers
% [st, re] = system('qstat -g c | grep matlab.q');
% [Tok, Rem] = strtok(re);
% [Tok, Rem] = strtok(Rem);
% [Tok, Rem] = strtok(Rem);
% [Tok, Rem] = strtok(Rem);
% [available] = strtok(Rem)
% parpool('ips_base', str2num(available))
%%%%% END CLUSTER PARPOOL  %%%%%%


%%%%%%%% EDIT %%%%%%%%%%%
GLMs =       {'block', ...
              'event'};
glmAnDirs =  {'analysis_block_acpc_lhPPC', ...
              'analysis_event_acpc_lhPPC'} 
% glmAnDirs =  {'analysis_block_acpc_lhIFG', ...
%   'analysis_event_acpc_lhIFG'} 
% glmAnDirs =  {'analysis_block_acpc_lhPPC', ...
%   'analysis_event_acpc_lhPPC'} 
conNums   =  {34, ...
              33}
Tmin = 0;  % Esto es igual a 0.001, 2.04 = 0.02, 1.6 = 0.05
% ultimo VOT = V15
% ultimo IFG = V01
% ultimo PPC = V02
versionNum = 'V02'
anArea = 'PPC'

% V04 =  testeo, borre todo lo anterior con referencias a SPM
% V05 = Genero excel V05 bien, pero no tienen info sobre fsaverage (305)
% V06 = Esta version con vtx sobre 305 pasa a analisis en R
% V07 = la he usado de pruebas junto con la 6 para el desarrollo de la
%       funcion myCreateProbOverlay.m. 
% V08 = habia bug, no buscaba bien dentro de roi pq vertex movido x 1.
%        Estos datos generan overlays en fsaverage y datos para excel. 
% V09: sera la ultima con todos los datos, a partir de ahora usare solo
%      unos contrastes selectos, ademas:    
%      -OK- hacer que cuando no haya T significativa guarde NA en vez de 0 en
%      toda la linea
%      -OK- ademas, cuando hago los mapas probabilisticos, hacerlo
%        thresholdeado
%      - OK- ademas usar el maximo y escribirlo en el filename para luego poder
%      pasarselo como parametro en el overlay (en myCreateProbOverlay.m)
%      -OK- hacer que siempre lea Contrasts.m desde la raiz, en fMRI_SPM
%      -OK- arreglar            el Bug que tiene yMin a la hora de hacer
%      probabilisticos, estaba usando el mismo surf que en el
%      esoapacio ind.
% V10: Despues de hacer tractografia, crear nuevas ROI-s
%      - Ahora voy a tener que cubrir todo el cortex con izda-dcha y
%         anterior-posterior. Que tamanos? Podria estar bien que el espacio
%         de busqueda sea el mismo o parecido que las termination points de
%         los tractos en Weiner-Yeatman. 

% 
%%%%%%%% END EDIT %%%%%%%

                
% Lo de a,c,pVWFA esta en myCreateLabels.m                
% VWFAstring = {'aVWFA', 'cVWFA', 'pVWFA'};
% dilateLabelBy = {'4', '8', '16', '32'};  % Ademas existe el 1 con un vertex
dilateLabelBy = {'16'};
% OJO CON LA CHAPUZA EN LA LINEA 182, ESTO LE AFECTA
if strcmp(anArea, 'PPC')
    VWFAletter = { 'lhNotVot_vof', 'lhNotVot_parc','lhPPC'};
elseif strcmp(anArea, 'IFG')
    VWFAletter = { 'lhNotVot_vof', 'lhNotVot_parc','lhIFG'};
else
    VWFAletter = { 'vof', 'parc','lhVOT'};
end
% lhVotNoV1V2yMin16.label  % En vez de usar GM y GMyMin usaremos este label


                
for noglm=1:length(GLMs)
    %%%%%%%%%%%%%%%%%
    glm      = GLMs{noglm}
    glmAnDir = glmAnDirs{noglm};
    conN   = conNums{noglm};
    
    basedir = [fMRIdir fsp glm fsp glmAnDir fsp 'SUBJECTS'];
    cd(basedir)
    groupscurly; % En el archivo groupscurly en SUBJECT guardar los sujetos que queramos analizar
    % quitar el S067, no tiene tractos
    allUnique4groupAnalysis = todos([1:61, 63:end]);
    
    addpath(fMRIdir);
    [nombres, numeros] = Contrasts(noglm);
    ContrastNames  = nombres;
    ContrasNumbers = numeros;
    rmpath(fMRIdir);
    % ContrastNames = TranslateConNames{noglm}; Lo he convertido a una func
    %%%%%%%%%%%%%%%%%
    
    parfor conNum = 1:length(numeros) 
        conNum4str = ContrasNumbers{conNum};
        conName    = ContrastNames{conNum};
        % Initialize structure
        maximas = struct('name', {}); 
        for ii = 1:length(allUnique4groupAnalysis)    
            maximas{ii}.name = allUnique4groupAnalysis{ii};
        end
        

        % Calculate values per every subject for this Contrast
        for ns = 1:length(maximas)
            subname = maximas{ns}.name

            % Read data: 
            TalXFM = xfm_read([fsdir fsp subname fsp 'mri' fsp 'transforms' ...
                               fsp 'talairach.xfm']);
            T1 = MRIread([fsdir fsp subname fsp 'mri' fsp 'T1.mgz']);
            % [aparcVert,aparcLabel,aparcClrtble] = read_annotation([fsdir ...
            %                    fsp subname ...
            %                    fsp 'label' ...
            %                    fsp 'lh.aparc.annot']);
            % V1thLabel = read_label(subname, 'lh.V1.thresh');
            % highres = MRIread([fMRIdir fsp glm fsp 'data' fsp subname ...
            %                    fsp 'anat' ...
            %                    fsp 'highres.nii']);
            lhwhite = read_surf([fsdir fsp subname fsp 'surf' fsp 'lh.white']);
            % lhpial = read_surf([fsdir fsp subname fsp 'surf' fsp 'lh.pial']);                           
            % spmT = MRIread([basedir fsp subname fsp 'results' ...
            %                    fsp 'spmT_' conNum4str '.img']);
            spmTsurf = MRIread([basedir fsp subname fsp 'results' ...
                               fsp conName '.mgh']);
            spmTsurf305 = MRIread([basedir fsp subname fsp 'results' ...
                                fsp conName '305.mgh']);            
            % Buscar la yMin
            yMin  = inv(TalXFM) * yMin305; 
            yMinSurf = T1.tkrvox2ras * inv(T1.vox2ras) * yMin;
            
            % Vamos a guardar solo el vertex, la T y la coord 152
            % INITIALIZE
%             prefijo = 'GM_';
%             maximas{ns}.([prefijo 'vtx'])    = NaN;
%             maximas{ns}.([prefijo 'vtx305']) = NaN;
%             maximas{ns}.([prefijo 'x'])      = NaN;
%             maximas{ns}.([prefijo 'y'])      = NaN;
%             maximas{ns}.([prefijo 'z'])      = NaN;
%             maximas{ns}.([prefijo 'T'])      = NaN;
%             
%             prefijo = 'GMyMin_';
%             maximas{ns}.([prefijo 'vtx'])    = NaN;
%             maximas{ns}.([prefijo 'vtx305']) = NaN;
%             maximas{ns}.([prefijo 'x'])      = NaN;
%             maximas{ns}.([prefijo 'y'])      = NaN;
%             maximas{ns}.([prefijo 'z'])      = NaN;
%             maximas{ns}.([prefijo 'T'])      = NaN;
            
            for jj = 1:length(VWFAletter) 
                for kk=1:length(dilateLabelBy)
                    prefijo = [VWFAletter{jj} dilateLabelBy{kk} '_']
                    maximas{ns}.([prefijo 'vtx'])    = NaN;
                    maximas{ns}.([prefijo 'vtx305']) = NaN;
                    maximas{ns}.([prefijo 'x'])      = NaN;
                    maximas{ns}.([prefijo 'y'])      = NaN;
                    maximas{ns}.([prefijo 'z'])      = NaN;
                    maximas{ns}.([prefijo 'fsx'])      = NaN;
                    maximas{ns}.([prefijo 'fsy'])      = NaN;
                    maximas{ns}.([prefijo 'fsz'])      = NaN;
                    maximas{ns}.([prefijo 'T'])      = NaN;
                    maximas{ns}.([prefijo 'inParc'])      = NaN;
                    maximas{ns}.([prefijo 'inVof'])      = NaN;
                end
            end
            
    
            %% GM sin yMin
            % ScannerRAS = T1.vox2ras*inv(T1.tkrvox2ras)*SurfaceRAS;
%             [spmTsurf_Max,spmTsurf_ind] = max(spmTsurf.vol(:));
%             [spmTsurf_Max305,spmTsurf_ind305] = max(spmTsurf305.vol(:));
%             if spmTsurf_Max >= 0.1 &&  spmTsurf_Max305 >= 0.1
%                 % Obtener la coordenada en la surface de lh.white:
%                 GM_spmTsurf_Surf = lhwhite(spmTsurf_ind, :);
%                 GMmghSurf        = T1.vox2ras*inv(T1.tkrvox2ras)*[GM_spmTsurf_Surf';1];
%                 % Al espacio MNI152
%                 GM_SPM305 = TalXFM      * GMmghSurf;
%                 GMmgh152  = MNI305to152 * GM_SPM305;
%                 % Write values
%                 prefijo = 'GM_';
%                 maximas{ns}.([prefijo 'vtx'])    = spmTsurf_ind - 1;
%                 maximas{ns}.([prefijo 'vtx305']) = spmTsurf_ind305 - 1;
%                 maximas{ns}.([prefijo 'x'])      = GMmgh152(1);
%                 maximas{ns}.([prefijo 'y'])      = GMmgh152(2);
%                 maximas{ns}.([prefijo 'z'])      = GMmgh152(3);
%                 maximas{ns}.([prefijo 'T'])      = spmTsurf_Max;
%             end
            
            %% GM con yMin
            % Primero hacemos copia y luego thresholdeamos
%             yMinspmTsurf = spmTsurf;
%             yMinspmTsurf.vol(lhwhite(:,2)  > yMinSurf(2)) = 0; 
%             yMinspmTsurf305 = spmTsurf305;
%             yMinspmTsurf305.vol(lhwhite305(:,2)  > yMin305(2)) = 0; 
%             % Obtain max and vertex
%             [yMinspmTsurf_Max, yMinspmTsurf_ind] = max(yMinspmTsurf.vol(:));
%             [yMinspmTsurf_Max305, yMinspmTsurf_ind305] = max(yMinspmTsurf305.vol(:));
%             if yMinspmTsurf_Max >= 0.1 && yMinspmTsurf_Max305 >= 0.1
%                 % Obtener la coordenada en la surface de lh.white:
%                 yMinGM_spmTsurf_Surf = lhwhite(yMinspmTsurf_ind, :);
%                 yMinGMmghSurf = T1.vox2ras*inv(T1.tkrvox2ras)*[yMinGM_spmTsurf_Surf';1];
%                 % En espacio MNI152
%                 yMinGM_SPM305 = TalXFM      * yMinGMmghSurf;
%                 yMinGMmgh152  = MNI305to152 * yMinGM_SPM305;
%                 % Write values
%                 prefijo = 'GMyMin_';
%                 maximas{ns}.([prefijo 'vtx'])    = yMinspmTsurf_ind - 1;
%                 maximas{ns}.([prefijo 'vtx305']) = yMinspmTsurf_ind305 - 1;
%                 maximas{ns}.([prefijo 'x'])      = yMinGMmgh152(1);
%                 maximas{ns}.([prefijo 'y'])      = yMinGMmgh152(2);
%                 maximas{ns}.([prefijo 'z'])      = yMinGMmgh152(3);
%                 maximas{ns}.([prefijo 'T'])      = yMinspmTsurf_Max; 
%             end
            
            %% LM-s dentro de los ROIs
            for jj = 1:length(VWFAletter) 
                for kk=1:length(dilateLabelBy)
                    % Read the ROI
                    % roiname = [VWFAletter{jj} 'VWFA' dilateLabelBy{kk}];
                    if strcmp([VWFAletter{jj} dilateLabelBy{kk}], 'lhVOT16')
                        roiname = 'lhVotNoV1V2yMin16';
                    else
                        roiname = [VWFAletter{jj} dilateLabelBy{kk}];
                    end
                    setenv('SUBJECTS_DIR', fsdir);
                    ROI    = read_label(subname, roiname);
                    % ojo, los roi-s creados por mi en el espacio fsaverage
                    % son el mismo para cada sujeto, pero el creado a
                    % partir de las manchas de DTI (vof o pArc) son
                    % diferentes para cada sujeto en el espacio de cada
                    % sujeto y en fsaverage
                    if nnz(ismember([1,2], jj))
                        ROI305 = read_label(subname, [roiname '_305']);
                    else
                        ROI305 = read_label('fsaverage', roiname);
                    end
                    % Copia y thresholdeo
                    tmpSurf = spmTsurf.vol;
                    tmpSurf(   1, setdiff(1:size(tmpSurf,   2), [(ROI(   :,1)+1)'])) = 0; % index es base 1 en matlab pero leo vertex en base 0 
                    tmpSurf305 = spmTsurf305.vol;
                    tmpSurf305(1, setdiff(1:size(tmpSurf305,2), [(ROI305(:,1)+1)'])) = 0; 
                    % Obtain max and vertex
                    [tmpSurf_Max, tmpSurf_ind] = max(tmpSurf(:));
                    [tmpSurf_Max305, tmpSurf_ind305] = max(tmpSurf305(:));
                    % Si no hay ningun voxel activo en el roi, dejar todo
                    % como NA, incluso en fsaverage, no siempre concuerda a
                    % la perfeccion, los dos tienen que ser no ceros
                    if tmpSurf_Max >= 0.1 && tmpSurf_Max305 >= 0.1
                        % Obtener la coordenada en la surface de lh.white:
                        tmpSurf_Surf = lhwhite(tmpSurf_ind, :);
                        RAS = T1.vox2ras*inv(T1.tkrvox2ras)*[tmpSurf_Surf';1];
                        % En espacio MNI152
                        RAS305 = TalXFM       * RAS;
                        RAS152  = MNI305to152 * RAS305;  
                        % Lo mismo partiendo de 305 (SurfaceRAS y RAS es =
                        % en fsaverage)
                        fsRAS305 = lhwhite305(tmpSurf_ind305, :);
                        fsRAS152  = MNI305to152 * [fsRAS305,1]'; 
                        % Vemos si la vtx305 esta dentro de vof y/o parc
                        parc = read_label(subname, 'parc16_305');
                        vof  = read_label(subname, 'vof16_305');
                        vofFlag = 0;
                        parcFlag = 0;
                        if (nnz(vof(:,1) == (tmpSurf_ind305 - 1))) > 0
                            vofFlag = 1;
                        end
                        if (nnz(parc(:,1) == (tmpSurf_ind305 - 1)))  > 0
                            parcFlag = 1;
                        end
                        % Write results
                        prefijo = [VWFAletter{jj} dilateLabelBy{kk} '_'];
                        maximas{ns}.([prefijo 'vtx'])    = tmpSurf_ind - 1;
                        maximas{ns}.([prefijo 'vtx305']) = tmpSurf_ind305 - 1;
                        maximas{ns}.([prefijo 'x'])      = RAS152(1);
                        maximas{ns}.([prefijo 'y'])      = RAS152(2);
                        maximas{ns}.([prefijo 'z'])      = RAS152(3);
                        maximas{ns}.([prefijo 'fsx'])      = fsRAS152(1);
                        maximas{ns}.([prefijo 'fsy'])      = fsRAS152(2);
                        maximas{ns}.([prefijo 'fsz'])      = fsRAS152(3);
                        maximas{ns}.([prefijo 'T'])      = tmpSurf_Max;
                        maximas{ns}.([prefijo 'inParc'])      = parcFlag;
                        maximas{ns}.([prefijo 'inVof'])      = vofFlag;
                    end
                end
            end            
            
            
 
    
        end
        
        
        
        % For every contrast write an excel with data of all subjects
        struct2csv(...
                cell2mat(maximas), ...
                char([basedir fsp glm 'Individual_acpc_' anArea '_con-' ...
                 conNum4str '_Tmin' num2str(Tmin) ...
                 '_yMin-40_' versionNum '.csv'])); 
        
        
        % Pintar los ROIs
        % myCreateProbOverlay(maximas, glm, conName, fsdir, basedir, versionNum, ...
        %                  Tmin, tmplate305, dilateBy)
%         myCreateProbOverlay(maximas, glm, conName, fsdir, basedir, versionNum, ...
%                             Tmin, spmTsurf305, '4')
        
        

            
    end

end

% subname = 'S006'
% [status, results] = system(['freeview -viewport 3d ' ...
%         '-f ' fsdir fsp subname fsp 'surf' fsp ...
%         'lh.inflated:annot=aparc.annot' ...
%         ':label=' fsdir fsp subname fsp 'label' fsp 'aVWFA16.label' ...
%         ':label=' fsdir fsp subname fsp 'label' fsp 'cVWFA16.label' ...
%         ':label=' fsdir fsp subname fsp 'label' fsp 'pVWFA16.label &' ...
%        ]);  





