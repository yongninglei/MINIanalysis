function [H,P,CI,STATS] = mySavefsavgSurfaceValues(trt, subs, TH)
% myqMRIAnalysisInROIs('TEST', TESTsubs, tempmgh)
    % trt = 'TEST'
    % subs = TESTsubs
    % trt = 'RETEST'
    % subs = DAY2subs


    fMRIconts = {'VOT_block_RWvsCB','VOT_block_RWvsCS','VOT_block_RWvsFF',...
                 'VOT_block_RWvsPS','VOT_block_RWvsPW','VOT_block_RWvsSD','VOT_block_RWvsNull'};
    for yy=1:length(fMRIconts)
        y = [subs.(fMRIconts{yy})]';
        datos = struct2cell(subs);
        datos = table(datos(1,:)');
        datos.Properties.VariableNames = {'SUBJECT'};
        datos = [datos, array2table(y)];
        filename = fullfile(MINIPath,'DATA','nogit_spmTvalues', [trt '_' fMRIconts{yy} '.csv']);
        % The tables are too slow to be read in R...
        % writetable(datos, ...
        %        filename, ...
        %        'FileType', 'text', ...
        %        'Delimiter', 'comma', ...
        %        'WriteVariableNames', true)  
        % In order to make the files smaller, create the mean of subjects here
        y = y';
        y(y<TH) = NaN;
        ym = mean(y,2, 'omitnan');
        csvwrite(filename, ym);
    end


end