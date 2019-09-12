function rootPath = rightMINIPath()
% Determine path to root of the mrVista directory
%
%        rootPath = vistaRootPath;
%
% This function MUST reside in the directory at the base of the
% MINI directory structure 
%
% Copyright Stanford team, mrVista, 2018

rootPath = which('rightMINIPath');

rootPath = fileparts(rootPath);

return
