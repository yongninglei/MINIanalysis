function myCreateDWIfsaverage(trt, subs, tempmgh)
    
    vOF    = [subs.DWI_vOF]';
    pARC   = [subs.DWI_pARC]';
    Th2label = 0.75;
    vOF(vOF    >= Th2label) = 1; % Reading the result of vol2surf and surf2surf, this way we binarize it
    pARC(pARC  >= Th2label) = 1;
    vOF(vOF    < Th2label) = 0; 
    pARC(pARC  < Th2label) = 0;
    
    vOFallperc  = 100 * (sum(vOF ,1,'default','omitnan') / size(vOF,1));
    pARCallperc = 100 * (sum(pARC,1,'default','omitnan') / size(pARC,1));
    
    % Tract intersection, the whole thing
    Th = 15;
    vOFallpercTh = vOFallperc;
    vOFallpercTh(vOFallpercTh <= Th) = 0;
    pARCallpercTh = pARCallperc;
    pARCallpercTh(pARCallpercTh <= Th) = 0;
    % Create the intersection     
    intersec  = sum([vOFallpercTh; pARCallpercTh],1,'default','omitnan') / 2;
    intersec(vOFallpercTh  == 0) = 0;
    intersec(pARCallpercTh == 0) = 0;
    
    % Write files    
    comp = tempmgh;
    comp.vol =  intersec;
    MRIwrite(comp, fullfile(MINIPath, 'DATA', 'DWI', [trt '_intersection.mgh']));
    comp = tempmgh;
    comp.vol =  pARCallpercTh;
    MRIwrite(comp, fullfile(MINIPath, 'DATA', 'DWI', [trt '_pARC.mgh']));
    comp = tempmgh;
    comp.vol =  vOFallpercTh;
    MRIwrite(comp, fullfile(MINIPath, 'DATA', 'DWI', [trt '_vOF.mgh']));
    
    
    
end