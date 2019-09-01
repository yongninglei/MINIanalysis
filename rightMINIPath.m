function rootPath = MINIPath()
% Determine path to root of the mrVista directory
%
%        rootPath = vistaRootPath;
%
% This function MUST reside in the directory at the base of the
% MINI directory structure 
%
% Copyright Stanford team, mrVista, 2018

rootPath = which('MINIPath');

rootPath = fileparts(rootPath);

return
