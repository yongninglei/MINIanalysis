function myDWIfMRIcorticalEndings(trt, subs)
% trt = 'TEST';
% subs = TESTsubs;

% Comparar valores funcionales dentro de los cortical endings vOF y pAF
vOF    = [subs.DWI_vOF]';
pARC   = [subs.DWI_pARC]';
Th2label = 0.75;
vOF(vOF    >= Th2label) = 1; % Reading the result of vol2surf and surf2surf, this way we binarize it
pARC(pARC  >= Th2label) = 1;
vOF(vOF    < Th2label) = 0; 
pARC(pARC  < Th2label) = 0;
% Leemos la versin label. Mirar, creo que al crear label thresholdeaba en 1 y
% aqui lo he hecho en 0.75
vOF_label    = [subs.DWI_vOF_label]';
pARC_label   = [subs.DWI_pARC_label]';
vOF_notVotlabel    = [subs.DWI_vOF_notVotlabel]';
pARC_notVotlabel   = [subs.DWI_pARC_notVotlabel]';

% Ahora asegurarnos de que no tiene ni V1V2 ni yMin
% DWI(isnan(vOF_label)) = 0;
% DWI(isnan(vOF)) = 0;
% DWI(isnan(pARC_label)) = 0;
% DWI(isnan(pARC)) = 0;
vOFvOT = vOF .* vOF_label;
pARCvOT = pARC .* pARC_label;
pARC2vOT  = pARCvOT;
pARC2vOT(pARC2vOT>0) = 2;
DWIvOT    = vOFvOT + pARC2vOT; % 0:ninguno, 1:vOF, 2:pARC, 3:los dos
DWIvOT(isnan(DWIvOT)) = 0;

% Obtener ahora valores medios para spm
thresh = 0.01;  % -inf, 0.01, 1.65
CB    = [subs.VOT_block_RWvsCB]'; CB(CB<=thresh) = NaN;
PS    = [subs.VOT_block_RWvsPS]'; PS(PS<=thresh) = NaN;
SD    = [subs.VOT_block_RWvsSD]'; SD(SD<=thresh) = NaN;
CS    = [subs.VOT_block_RWvsCS]'; CS(CS<=thresh) = NaN;
FF    = [subs.VOT_block_RWvsFF]'; FF(FF<=thresh) = NaN;
PW    = [subs.VOT_block_RWvsPW]'; PW(PW<=thresh) = NaN;
% Initialize NaN matrix for the averages of PER and WL
PER   = NaN * ones(size(CB));
WL    = NaN * ones(size(CB));
% Initialize NaN vector that should have one value per subject depending on
% tract
PER_vOF = NaN * ones(size(CB,1),1);
PER_pAF = NaN * ones(size(CB,1),1);
PER_INT = NaN * ones(size(CB,1),1);
WL_vOF  = NaN * ones(size(CB,1),1);
WL_pAF  = NaN * ones(size(CB,1),1);
WL_INT  = NaN * ones(size(CB,1),1);
% First create the average PER or WL per subject,
% Convert the NaN-s to zeros (so that don't sum in the next matrix multiplication
% Then create the average of all values that are vOF, or pAF or INT for that
% subject, using a matrix multiplication / number of nonzero values
for nsub=1:size(CB,1)
    tmp = mean([CB(nsub,:);PS(nsub,:);SD(nsub,:)],1,'omitnan');
    tmp(isnan(tmp)) = 0;
    PER(nsub,:)     = tmp;
    
    PER_vOF(nsub) = (PER(nsub,:) * (DWIvOT(nsub,:)==1)')/nnz(PER(nsub,:).* (DWIvOT(nsub,:)==1));
    PER_pAF(nsub) = (PER(nsub,:) * (DWIvOT(nsub,:)==2)')/nnz(PER(nsub,:).* (DWIvOT(nsub,:)==2));
    PER_INT(nsub) = (PER(nsub,:) * (DWIvOT(nsub,:)==3)')/nnz(PER(nsub,:).* (DWIvOT(nsub,:)==3));
end
for nsub=1:size(CB,1)
    tmp = mean([CS(nsub,:);FF(nsub,:);PW(nsub,:)],1,'omitnan');
    tmp(isnan(tmp)) = 0;
    WL(nsub,:) = tmp;
    
    WL_vOF(nsub) = (WL(nsub,:) * (DWIvOT(nsub,:)==1)')/nnz(WL(nsub,:).* (DWIvOT(nsub,:)==1));
    WL_pAF(nsub) = (WL(nsub,:) * (DWIvOT(nsub,:)==2)')/nnz(WL(nsub,:).* (DWIvOT(nsub,:)==2));
    WL_INT(nsub) = (WL(nsub,:) * (DWIvOT(nsub,:)==3)')/nnz(WL(nsub,:).* (DWIvOT(nsub,:)==3));
end


% Remove the subjects we know there are bad
% Make this code in tables and remove these subjects
if strcmp(trt, 'TEST')
    removeSUBS   = ["S013","S018","S004","S029","S032","S048","S067"];
    removeSUBS   = [4,13,18,29,32,43,62];
else
    removeSUBS = ["S072", "S086", "S097"];
    2,16
end

doRemoveSUBS = true;
if doRemoveSUBS
    WL_pAF(removeSUBS)  = NaN;
    WL_vOF(removeSUBS)  = NaN;
    WL_INT(removeSUBS)  = NaN;
    PER_pAF(removeSUBS) = NaN;
    PER_INT(removeSUBS) = NaN;
    PER_vOF(removeSUBS) = NaN;
end

mean([WL_pAF,WL_INT,WL_vOF,PER_pAF,PER_INT,PER_vOF],1,'omitnan');
disp('LEX, pAF vs vOF')
[H,P,CI,STATS]= ttest(WL_pAF,WL_vOF,'tail','right')
disp('PER, pAF vs vOF')
[H,P,CI,STATS]= ttest(PER_pAF,PER_vOF,'tail','left')

% Meter ademas la intersection a ver si vemos escalado:
% % WL
% [H,P,CI,STATS]= ttest(WL_pAF,WL_INT,'tail','right')
% [H,P,CI,STATS]= ttest(WL_INT,WL_vOF,'tail','right')
% % PER
% [H,P,CI,STATS]= ttest(PER_pAF,PER_INT,'tail','left')
% [H,P,CI,STATS]= ttest(PER_INT,PER_vOF,'tail','left')
              