function myDWIfMRIDicotomicData(trt, subs, kkvertex)
    
  
    
    vOF    = [subs.DWI_vOF]';
    pAF    = [subs.DWI_pARC]';
    Th2label = 0.75;
    vOF(vOF    >= Th2label) = 1; % Reading the result of vol2surf and surf2surf, this way we binarize it
    pAF(pAF    >= Th2label) = 1;
    vOF(vOF    < Th2label)  = 0; 
    pAF(pAF    < Th2label)  = 0;
    pAF_aVWFA = pAF(:, kkvertex.aVWFA');
    pAF_pVWFA = pAF(:, kkvertex.pVWFA');
    vOF_aVWFA = vOF(:, kkvertex.aVWFA');
    vOF_pVWFA = vOF(:, kkvertex.pVWFA');

    pAF_mOTS = pAF(:, kkvertex.([trt '_lexVWFA4'])');
    pAF_pOTS = pAF(:, kkvertex.([trt '_perVWFA4'])');
    vOF_mOTS = vOF(:, kkvertex.([trt '_lexVWFA4'])');
    vOF_pOTS = vOF(:, kkvertex.([trt '_perVWFA4'])');

    % Create the dicotomic variable, first sum to be able to threshold afterwards
    sum_pAF_aVWFA = sum(pAF_aVWFA, 2, 'native','omitnan');
    sum_pAF_pVWFA = sum(pAF_pVWFA, 2, 'native','omitnan');
    sum_vOF_aVWFA = sum(vOF_aVWFA, 2, 'native','omitnan');
    sum_vOF_pVWFA = sum(vOF_pVWFA, 2, 'native','omitnan');

    sum_pAF_mOTS = sum(pAF_mOTS, 2, 'native','omitnan');
    sum_pAF_pOTS = sum(pAF_pOTS, 2, 'native','omitnan');
    sum_vOF_mOTS = sum(vOF_mOTS, 2, 'native','omitnan');
    sum_vOF_pOTS = sum(vOF_pOTS, 2, 'native','omitnan');

    % Create the dicotomic one
    dic_pAF_aVWFA = (sum_pAF_aVWFA > 0);
    dic_pAF_pVWFA = (sum_pAF_pVWFA > 0);
    dic_vOF_aVWFA = (sum_vOF_aVWFA > 0);
    dic_vOF_pVWFA = (sum_vOF_pVWFA > 0);

    dic_pAF_mOTS = (sum_pAF_mOTS > 0);
    dic_pAF_pOTS = (sum_pAF_pOTS > 0);
    dic_vOF_mOTS = (sum_vOF_mOTS > 0);
    dic_vOF_pOTS = (sum_vOF_pOTS > 0);

    % Create a table:
    mat2R_lit   = [dic_pAF_aVWFA,dic_pAF_pVWFA,dic_vOF_aVWFA,dic_vOF_pVWFA];
    mat2R_MINI  = [dic_pAF_mOTS,dic_pAF_pOTS,dic_vOF_mOTS,dic_vOF_pOTS];
    datos = array2table([mat2R_lit, mat2R_MINI]);
    datos.Properties.VariableNames = {'pAF_aVWFA','pAF_pVWFA','vOF_aVWFA','vOF_pVWFA', ...
                                      'pAF_mOTS' ,'pAF_pOTS' ,'vOF_mOTS' ,'vOF_pOTS'};
    tmpsubs = struct2cell(subs);
    datos.SUBJECT = tmpsubs(1,:)';

    % Write it    
    filename = [trt '_DWIintoROIs_litVWFA_and_mOTS-pOTS.csv'];
    writetable(datos, ...
               fullfile(MINIPath,'DATA','DWI',filename), ...
               'FileType', 'text', ...
               'Delimiter', 'comma', ...
               'WriteVariableNames', true)
    
    
end