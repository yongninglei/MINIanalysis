function findMissingFoldersAndFiles(basedir, subjectIds, contrasts, sm)
    
    % Initialize the file to store missing folders or files
    missingFilePath = fullfile(basedir, 'missingFoldersAndFiles.txt');
    % Open or create the file for writing
    fid = fopen(missingFilePath, 'w');
    
    if fid == -1
        error('Unable to create or open the missing file list text file.');
    end
    
    % Loop through each subject ID
    for i = 1:length(subjectIds)
        subjectFolderName=sprintf('S%03d', subjectIds(i));
        subjectFolder = fullfile(basedir, subjectFolderName);
        fprintf(fid, 'Searching under folder: %s\n', subjectFolder);
        if ~exist(subjectFolder, 'dir')
            % Subject folder does not exist
            fprintf(fid, 'Missing folder: %s\n', subjectFolderName);
            fprintf('Missing folder: %s\n', subjectFolderName);
        else
            % Check for each specified contrast file
            for j = 1:length(contrasts)
                filePatterns = {sprintf('%s305%s.mgh', contrasts{j}, sm), sprintf('%s305.mgh', contrasts{j})};
                for k = 1:length(filePatterns)
                    filePath = fullfile(subjectFolder, 'results', filePatterns{k});
                    if ~exist(filePath, 'file')
                        % File does not exist
                        fprintf(fid, 'Missing file: %s under %s\n', filePatterns{k}, subjectFolderName);
                        fprintf('Missing file: %s under %s\n', filePatterns{k}, subjectFolderName);
                    end
                end
            end
        end
    end
    
    % Close the file
    fclose(fid);
end
