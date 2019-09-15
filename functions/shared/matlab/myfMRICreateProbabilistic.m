function myfMRICreateProbabilistic(trt, subs, kkvertex,tempmgh,sm)
    
  contrastes   = {'RWvsCB','RWvsPS','RWvsSD','RWvsCS','RWvsFF','RWvsPW','RWvsNull'};
% ClustSize = 100;  % hacerlo dinamico, coge el 60% por ejemplo
prcntj = .5;
binbrain = {};

for contr=1:length(contrastes)
    todos    = [subs.(['VOT_block_' contrastes{contr}])]';
    todos    = todos(:, kkvertex.adjLOTS);
    brain    = repmat(tempmgh.vol, [size(todos,1),1]);
    brain(:, kkvertex.adjLOTS) = todos;
    
    % Lo tengo que hacer por sujeto independientemente
    for ii=1:size(brain,1)
        [sortedBrain,sortingIndices] = sort(brain(ii,:),'descend');
        Maximo = sortedBrain(1);
        percMaximo = prcntj * Maximo;
        % percMaximo = 1.65;
        ClustSize = size(find(sortedBrain > percMaximo),2);
        nonMaxValueIndices = sortingIndices(ClustSize+1:end);
        brain(ii,nonMaxValueIndices) = 0;
    end
    % Make binary sums
        binbrain{contr} =  (brain > 0);
    % Make Averages
        % paraEscribir = mean(brain,1,'omitnan');
        
    % WRITE MGH
    if contr == 3
       perceptual = sum([binbrain{1};binbrain{2};binbrain{3}], ...
                         1,'omitnan');
       paraEscribir = 100*(perceptual/(3*size(brain,1)));
       logP = tempmgh;
       logP.vol = paraEscribir;
       MRIwrite(logP, fullfile(MINIPath,'DATA','fMRI',[trt '_PERCEPTUAL_adjLOTS2_ProbMap_Prcj' num2str(prcntj) sm '.mgh']));  
    elseif contr == 6
       semantic = sum([binbrain{4};binbrain{5};binbrain{6}], ...
                         1,'omitnan');
       paraEscribir = 100*(semantic/(3*size(brain,1)));
       logP = tempmgh;
       logP.vol = paraEscribir;
       MRIwrite(logP, fullfile(MINIPath,'DATA','fMRI',[trt '_LEXICAL_adjLOTS2_ProbMap_Prcj' num2str(prcntj) sm '.mgh']));  
    elseif contr == 7
       semantic = sum(binbrain{7},1,'omitnan');
       paraEscribir = 100*(semantic/size(brain,1));
       logP = tempmgh;
       logP.vol = paraEscribir;
       MRIwrite(logP, fullfile(MINIPath,'DATA','fMRI',[trt '_NULL_adjLOTS2_ProbMap_Prcj' num2str(prcntj) sm '.mgh']));  
    end
end

    
    
end