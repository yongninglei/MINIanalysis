function [l] = myFSread_label(sname, lname, use_path)
% l = myFSread_label(<sname>, lname)
%
% reads the label file 'lname' from the subject 'sname' 
% in the subject's label directory into the vector l
% l will be nvertices-by-5, where each column means:
% (1) vertex number, (2-4) xyz at each vertex, (5) stat
%
% IMPORTANT: the vertex number is 0-based.
% 
% NOTE GLU: added use_path to make it read the whole path, not in /label


%
% read_label.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.7 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


l = [];

if(nargin ~= 3)
  fprintf('l = read_label(<sname>, lname, use_path)\n');
  return;
end

if use_path
    fname = lname;
else
    if(~isempty(sname))
      sdir = getenv('SUBJECTS_DIR') ;
      fname = sprintf('%s/%s/label/%s.label', sdir, sname, lname) ;
    else
        ind = findstr(lname, '.label') ;
        if (length(ind) > 0)
            fname = lname ;
        else
            fname = sprintf('%s.label', lname);
        end
    end
end
% open it as an ascii file
fid = fopen(fname, 'r') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

fgets(fid) ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

line = fgets(fid) ;
nv = sscanf(line, '%d') ;
l = fscanf(fid, '%d %f %f %f %f\n') ;
l = reshape(l, 5, nv) ;
l = l' ;

fclose(fid) ;

