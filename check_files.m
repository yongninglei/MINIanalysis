% this is for cheking icloud data is missing or not
% specify the zip dir
zip_dir='/bcbl/home/public/Gari/MINI/DATA/icloud_data/MINI/all_the_zips_from_gDrive';
zip_files=dir(fullfile(zip_dir,'*zip'));

% create a temorary folder that check all the files
temp_dir=fullfile(zip_dir, 'temp_extraction');
if ~exist(temp_dir, "dir")
    mkdir(temp_dir);
end

% loop through all the zip files:

for i = i:length(zip_files)
    zip_file_path=fullfile(zip_files(i).folder, zip_files(i).name);
    fprintf('Checking the zip file %s\n', zip_files(i).name)
    folder= strrep(zip_file_path, '.zip', '');
    
    bolck_dir_path= fullfile(folder, 'MINI/DATA/myGLMFIT/fMRI_VOT/block/');
    
end
