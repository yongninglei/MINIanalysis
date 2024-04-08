 function subs = myReadStructures(subs, DATAdir, sm, tempmgh, SHOW, hemi)
fsp = filesep;
if strcmp(hemi,'rh')
    H = 'R';
else
    H = 'L';
end
for ns = 1: length(subs)
    sub = subs(ns).name;
    fnames = fieldnames(subs(ns));
    
    for ii = 2:length(fnames)
        fname = fnames{ii};
        switch fname
            case {'CT'}
                if SHOW disp('CT: reading the data... '); end;
                %% there is a typo in sm, sm is fhwm5, but in fsfoler, it fwhm5
                filename = [DATAdir fsp 'freesurferacpc' fsp sub fsp 'surf' fsp ...
                                hemi '.thickness.fwhm5.fsaverage.mgh'];
                if exist(filename)
                    temp = MRIread(filename);
                    subs(ns).(fname) = temp.vol';
                    if SHOW disp(' ...done');end;
                else
                    disp([sub ' ' fname ': file doesnt exist'])
                end
                
            case{'CURV'}
                if SHOW disp('CURV: reading data... ');end;
                filename  = [DATAdir fsp 'freesurferacpc' fsp sub fsp 'surf' fsp ...
                                hemi '.curv.fwhm5.fsaverage.mgh'];
                if exist(filename)
                    temp = MRIread(filename);
                    subs(ns).(fname) = temp.vol';
                    if SHOW disp(' ...done');end;
                else
                    disp([sub ' ' fname ': file doesnt exist'])
                end

%             case{'DWI_vOF', 'DWI_pARC'}
%                 if SHOW disp(['DWI: reading ' tract '... ']); end;
%                 tract = upper(fname(5:end));
%                 if strcmp(tract, 'VOF')
%                     tname = [H '_VOF_tracts305.mgh'];
%                 elseif strcmp(tract, 'PARC')
%                     tname = [H '_Arcuate_Posterior_tracts305.mgh'];
%                 end
%                 filename = [DATAdir fsp fname fsp sub fsp 'dmri' fsp tname];
%                     
%                 if exist(filename)
%                     temp = MRIread(filename);
%                     subs(ns).(fname) = temp.vol';
%                     if SHOW disp(' ...done');end;
%                 else
%                     disp([sub ' ' fname ': file doesnt exist'])
%                 end
%             case{'DWI_vOF_label',       'DWI_pARC_label', ...
%                  'DWI_vOF_notVotlabel', 'DWI_pARC_notVotlabel'}
%                 if SHOW disp(['DWI: reading ' tract '... ']); end;
%                 if H=='L'
%                     if contains(fname, 'vOF_label')
%                         tract = 'vof';
%                         pre   = 'NoV1V2_';
%                     elseif contains(fname, 'vOF_notVotlabel')
%                         tract = 'vof';
%                         pre   = 'lhNotVot_';
%                     elseif contains(fname, 'pARC_label')
%                         tract = 'parc';
%                         pre   = 'NoV1V2_';
%                     elseif contains(fname, 'pARC_notVotlabel')
%                         tract = 'parc';
%                         pre   = 'lhNotVot_';
%                     end
%                 else
%                     if contains(fname, [hemi '_vOF_label'])
%                         tract = [hemi '_vof'];
%                         pre   = [hemi '_NoV1V2_'];
%                     elseif contains(fname, [hemi '_vOF_notVotlabel'])
%                         tract = [hemi '_vof'];
%                         pre   = [hemi '_NotVot_'];
%                     elseif contains(fname, ['pARC_label'])
%                         tract = [hemi '_parc'];
%                         pre   = [hemi '_NoV1V2_'];
%                     elseif contains(fname, [hemi '_pARC_notVotlabel'])
%                         tract = [hemi '_parc'];
%                         pre   = [hemi '_NotVot_'];
%                     end
%                 end
%                 filename = [DATAdir fsp 'DWI' fsp sub fsp ...
%                                       'label' fsp pre tract '16_305.label'];
%                 if exist(filename,'file')
%                     temp = myFSread_label(sub, filename,1);
%                     tempall = tempmgh;
%                     tempall.vol = zeros(size(tempall.vol));
%                     tempall.vol((temp(:,1)+1)) = 1;
%                     subs(ns).(fname) = tempall.vol';
%                     if SHOW, disp(' ...done');end;
%                 else
%                     disp([sub ' ' fname ': file doesnt exist'])
%                 end
                                
            case {'VOT_block_RWvsCB','VOT_block_RWvsCS','VOT_block_RWvsFF', ...
                  'VOT_block_RWvsPS','VOT_block_RWvsPW','VOT_block_RWvsSD', ...
                  'VOT_event_RWvsCB','VOT_event_RWvsCS','VOT_event_RWvsFF', ...
                  'VOT_event_RWvsPS','VOT_event_RWvsPW','VOT_event_RWvsSD', ...
                  'VOT_event_RWvsNull','VOT_block_RWvsNull'}
                if SHOW, disp([fname '... ']);end
                design = fname(5:9);
                cont = fname(11:end);
                if strcmp(hemi,'lh')
                    %filename = [DATAdir fsp hemi '_fMRI_VOT' fsp design fsp sub ...
                    %                   fsp 'results' fsp cont '305' sm '.mgh'];
                    % change this only for tiger's testing using 
                    filename = ['/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/block/analysis_block_acpc_lhITfusLatOcc/smoothed_fhmw5'...
                                fsp sub fsp 'results' fsp cont '305' sm '.mgh'];
                else 
                    filename = [DATAdir fsp 'fMRI_VOT' fsp design fsp sub ...
                                       fsp 'results' fsp cont '305' sm '.mgh'];
                end
                
                if exist(filename,'file')
                    temp = MRIread(filename);
                    subs(ns).(fname) = temp.vol';
                    if SHOW, disp(' ...done');end
                else
                    disp([sub ' ' fname ': file doesnt exist'])
                end
%             case {'PPC_block_RWvsCB','PPC_block_RWvsCS','PPC_block_RWvsFF', ...
%                   'PPC_block_RWvsPS','PPC_block_RWvsPW','PPC_block_RWvsSD', ...
%                   'PPC_event_RWvsCB','PPC_event_RWvsCS','PPC_event_RWvsFF', ...
%                   'PPC_event_RWvsPS','PPC_event_RWvsPW','PPC_event_RWvsSD', ...
%                   'PPC_event_RWvsNull','PPC_block_RWvsNull'}
%                 if SHOW disp([fname '... ']);end;
%                 design = fname(5:9);
%                 cont = fname(11:end);
%                 filename = [DATAdir fsp 'fMRI_PPC' fsp design fsp sub ...
%                                        fsp 'results' fsp cont '305' sm '.mgh'];
%                 if exist(filename)
%                     temp = MRIread(filename);
%                     subs(ns).(fname) = temp.vol';
%                     if SHOW disp(' ...done');end;
%                 else
%                     disp([sub ' ' fname ': file doesnt exist'])
%                 end
%             case {'IFG_block_RWvsCB','IFG_block_RWvsCS','IFG_block_RWvsFF', ...
%                   'IFG_block_RWvsPS','IFG_block_RWvsPW','IFG_block_RWvsSD', ...
%                   'IFG_event_RWvsCB','IFG_event_RWvsCS','IFG_event_RWvsFF', ...
%                   'IFG_event_RWvsPS','IFG_event_RWvsPW','IFG_event_RWvsSD', ...
%                   'IFG_event_RWvsNull','IFG_block_RWvsNull'}
%                 if SHOW disp([fname '... ']);end;
%                 design = fname(5:9);
%                 cont = fname(11:end);
%                 filename = [DATAdir fsp 'fMRI_IFG' fsp design fsp sub ...
%                                        fsp 'results' fsp cont '305' sm '.mgh'];
%                 if exist(filename)
%                     temp = MRIread(filename);
%                     subs(ns).(fname) = temp.vol';
%                     if SHOW disp(' ...done');end;
%                 else
%                     disp([sub ' ' fname ': file doesnt exist'])
%                 end
% 
            case{'qMRI_MTV', 'qMRI_T1qMRI'}
                qM = fname(6:end);
                if SHOW, disp(['qMRI: reading ' qM '... ']);end;
                filename = ['/bcbl/home/public/Gari/MINI/DATA/icloud_data' fsp 'qMRI' fsp sub fsp ...
                            'OutPutFiles_1' fsp 'BrainMaps' fsp 'lh.' qM '_305' sm '.mgh'];
                if exist(filename,'file')
                    temp = MRIread(filename);
                    subs(ns).(fname) = temp.vol';
                    if SHOW disp(' ...done');end;
                else
                    disp([sub ' ' fname ': file doesnt exist'])
                end
%                 
%             case{'qMRI_MTV_WM', 'qMRI_T1qMRI_WM'}
%                 qM = fname(6:end);
%                 if SHOW, disp(['qMRI_WM: reading ' qM '... ']);end;
%                 filename = [DATAdir fsp 'qMRI_WM' fsp sub fsp ...
%                             'OutPutFiles_1' fsp 'BrainMaps' fsp qM '_305' sm '.mgh'];
%                 if exist(filename,'file')
%                     temp = MRIread(filename);
%                     subs(ns).(fname) = temp.vol';
%                     if SHOW disp(' ...done');end;
%                 else
%                     disp([sub ' ' fname ': file doesnt exist'])
%                 end
%                 
%             case {'VOL'}
%                 if SHOW, disp('VOL: reading the data... ');end;
%                 filename = [DATAdir fsp fname fsp sub fsp 'surf' fsp ...
%                                 'lh.volume.fsaverage' sm '.mgh'];
%                 if exist(filename,'file')
%                     temp = MRIread(filename);
%                     subs(ns).(fname) = temp.vol';
%                     if SHOW disp(' ...done');end;
%                 else
%                     disp([sub ' ' fname ': file doesnt exist'])
%                 end
            otherwise
                disp([fname ': Unknown field'])
        end    
    end
end
end