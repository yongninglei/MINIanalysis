function subs = myInitStructures(subs, tempmgh, fMRIareas, designs, Contrasts)
    % todos los datos con vectores con NaN, asi si luedo no existe este
% dato para este sujeto que no haga nada y saldra de los analisis posteriores. 
NaNVector = NaN(size(tempmgh.vol'));
for ns = 1: length(subs)
    sub = subs(ns).name;
    %subs(ns).CT                  = NaNVector;
    %subs(ns).CURV                = NaNVector;
    %subs(ns).DWI_vOF             = NaNVector;
    %subs(ns).DWI_pARC            = NaNVector;
    %subs(ns).DWI_vOF_label       = NaNVector;
    %subs(ns).DWI_pARC_label      = NaNVector;
    %subs(ns).DWI_vOF_notVotlabel = NaNVector;
    %subs(ns).DWI_pARC_notVotlabel= NaNVector;
    %subs(ns).qMRI_MTV            = NaNVector;
    %subs(ns).qMRI_T1qMRI         = NaNVector;
    %subs(ns).qMRI_MTV_WM         = NaNVector;
    %subs(ns).qMRI_T1qMRI_WM      = NaNVector;
    for area = fMRIareas; for design = designs; for contrast = Contrasts
        fname = [area{:} '_' design{:} '_' contrast{:}];
        subs(ns).(fname)    = NaNVector;
    end; end; end
    %subs(ns).VOL                 = NaNVector;
end
end