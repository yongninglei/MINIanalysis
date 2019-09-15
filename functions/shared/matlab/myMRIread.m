% This is my own version of MRIread from Freesurfer.
% I made the following modifications:
% 1. It doesn't crash if MRI parameters (TE/TR) are not available, just
% throws a warning instead.
% 2. It is parfor-friendly
% 3. It uses /cluster/scratch/monday as temporary directory (/tmp/ runs out
% of space sometimes...). This can be easily changed, though :-)
% Juan Eugenio Iglesias
% myMRIread(filename,justheader)
function MRI = myMRIread(filename,justheader)

tempdir='/scratch/';
if exist(tempdir,'dir')==0
    tempdir='/tmp/';
end
if exist(tempdir,'dir')==0
    tempdir='/cluster/scratch/monday/';
end

if ~isdeployed
    addpath('/usr/local/freesurfer/stable5_1_0/matlab');
    addpath('/autofs/homes/002/iglesias/matlab/myFunctions/');
end

if(exist('justheader')~=1) justheader = 0; end
if length(filename)>6 && strcmp(filename(end-6:end),'.nii.gz')>0
    aux=clock();
    aux2=num2str(round(1e6*aux(end)));
    aux=filename;
    aux(aux=='/')='_';
    aux(aux=='.')='_';
    fname2 = ['/tmp/' getenv('PBS_JOBID') '_' aux '_' aux2 '.nii'];
    while exist(fname2,'file')>0
        tic
        pause(rand(1));
        a=toc;
        fname2 = ['/tmp/' getenv('PBS_JOBID') '_' aux '_' num2str(round(a*1000000)) '_' aux2 '.mgh'];
    end
    system(['gunzip -c ' filename ' > ' fname2]);
    MRI=myMRIread_aux(fname2,justheader);
    delete(fname2);
elseif length(filename)>3 && strcmp(filename(end-3:end),'.mgz')>0
    aux=clock();
    aux2=num2str(round(1e6*aux(end)));
    aux=filename;
    aux(aux=='/')='_';
    aux(aux=='.')='_';
    fname2 = ['/tmp/' getenv('PBS_JOBID') '_' aux '_' aux2 '.mgh'];
    while exist(fname2,'file')>0
        tic
        pause(rand(1));
        a=toc;
        fname2 = ['/tmp/' getenv('PBS_JOBID') '_' aux '_' num2str(round(a*1000000)) '_' aux2 '.mgh'];
    end
    system(['gunzip -c ' filename ' > ' fname2]);
    MRI=myMRIread_aux(fname2,justheader);
    delete(fname2);
else
    MRI=myMRIread_aux(filename,justheader);
end





function mri = myMRIread_aux(fstring,justheader)


mri = [];

[fspec fstem fmt] = MRIfspec(fstring);
if(isempty(fspec))
    err = sprintf('ERROR: cannot determine format of %s (%s)\n',fstring,mfilename);
    error(err);
    return;
end

mri.srcbext = '';    % empty be default
mri.analyzehdr = []; % empty be default
mri.bhdr = []; % empty be default

%-------------- MGH ------------------------%
switch(fmt)
    case {'mgh'}
        [mri.vol, M, mr_parms, volsz] = load_mgh(fspec,[],[],justheader);
        if(isempty(M))
            fprintf('ERROR: loading %s as MGH\n',fspec);
            mri = [];
            return;
        end
        if(~justheader)
            mri.vol = permute(mri.vol,[2 1 3 4]);
            volsz = size(mri.vol);
        else
            mri.vol = [];
            volsz(1:2) = [volsz(2) volsz(1)];
        end
        try
            tr = mr_parms(1);
            flip_angle = mr_parms(2);
            te = mr_parms(3);
            ti = mr_parms(4);
        catch
            tr=nan; flip_angle=nan; te=nan; ti=nan;
        end
        %--------------- bshort/bfloat -----------------------%
    case {'bhdr'}
        if(~justheader)
            [mri.vol bmri] = fast_ldbslice(fstem);
            if(isempty(mri.vol))
                fprintf('ERROR: loading %s as bvolume\n',fspec);
                mri = [];
                return;
            end
            volsz = size(mri.vol);
        else
            mri.vol = [];
            bmri = fast_ldbhdr(fstem);
            if(isempty(bmri))
                fprintf('ERROR: loading %s as bvolume\n',fspec);
                mri = [];
                return;
            end
            [nslices nrows ncols ntp] = fmri_bvoldim(fstem);
            volsz = [nrows ncols nslices ntp];
        end
        [nrows ncols ntp fs ns endian bext] = fmri_bfiledim(fstem);
        mri.srcbext = bext;
        M = bmri.T;
        tr = bmri.tr;
        flip_angle = bmri.flip_angle;
        te = bmri.te;
        ti = bmri.ti;
        mri.bhdr = bmri;
        %------- analyze -------------------------------------
    case {'img'}
        hdr = load_analyze(fspec,justheader);
        if(isempty(hdr))
            fprintf('ERROR: loading %s as analyze\n',fspec);
            mri = [];
            return;
        end
        volsz = hdr.dime.dim(2:end);
        indnz = find(volsz~=0);
        volsz = volsz(indnz);
        volsz = volsz(:)'; % just make sure it's a row vect
        if(~justheader) mri.vol = permute(hdr.vol,[2 1 3 4]);
        else            mri.vol = [];
        end
        volsz([1 2]) = volsz([2 1]); % Make consistent. No effect when rows=cols
        tr = 1000*hdr.dime.pixdim(5); % msec
        flip_angle = 0;
        te = 0;
        ti = 0;
        hdr.vol = []; % already have it above, so clear it
        M = vox2ras_1to0(hdr.vox2ras);
        mri.analyzehdr = hdr;
        %------- nifti nii -------------------------------------
    case {'nii'}
        hdr = load_nifti(fspec,justheader);
        if(isempty(hdr))
            fprintf('ERROR: loading %s as analyze\n',fspec);
            mri = [];
            return;
        end
        volsz = hdr.dim(2:end);
        indnz = find(volsz~=0);
        volsz = volsz(indnz);
        volsz = volsz(:)'; % just make sure it's a row vect
        if(~justheader) mri.vol = permute(hdr.vol,[2 1 3 4]);
        else            mri.vol = [];
        end
        volsz([1 2]) = volsz([2 1]); % Make consistent. No effect when rows=cols
        tr = hdr.pixdim(5); % already msec
        flip_angle = 0;
        te = 0;
        ti = 0;
        hdr.vol = []; % already have it above, so clear it
        M = hdr.vox2ras;
        mri.niftihdr = hdr;
        %---------------------------------------------------
    otherwise
        fprintf('ERROR: format %s not supported\n',fmt);
        mri = [];
        return;
end
%--------------------------------------%

mri.fspec = fspec;
mri.pwd = pwd;

mri.flip_angle = flip_angle;
mri.tr  = tr;
mri.te  = te;
mri.ti  = ti;

% Assumes indices are 0-based. See vox2ras1 below for 1-based.  Note:
% MRIwrite() derives all geometry information (ie, direction cosines,
% voxel resolution, and P0 from vox2ras0. If you change other geometry
% elements of the structure, it will not be reflected in the output
% volume. Also note that vox2ras still maps Col-Row-Slice and not
% Row-Col-Slice.  Make sure to take this into account when indexing
% into matlab volumes (which are RCS).
mri.vox2ras0 = M;

% Dimensions not redundant when using header only
volsz(length(volsz)+1:4) = 1; % Make sure all dims are represented
mri.volsize = volsz(1:3); % only spatial components
mri.height  = volsz(1);   % Note that height (rows) is the first dimension
mri.width   = volsz(2);   % Note that width (cols) is the second dimension
mri.depth   = volsz(3);
mri.nframes = volsz(4);

%--------------------------------------------------------------------%
% Everything below is redundant in that they can be derivied from
% stuff above, but they are in the MRI struct defined in mri.h, so I
% thought I would add them here for completeness.  Note: MRIwrite()
% derives all geometry information (ie, direction cosines, voxel
% resolution, and P0 from vox2ras0. If you change other geometry
% elements below, it will not be reflected in the output volume.

mri.vox2ras = mri.vox2ras0;

mri.nvoxels = mri.height * mri.width * mri.depth; % number of spatial voxles
mri.xsize = sqrt(sum(M(:,1).^2)); % Col
mri.ysize = sqrt(sum(M(:,2).^2)); % Row
mri.zsize = sqrt(sum(M(:,3).^2)); % Slice

mri.x_r = M(1,1)/mri.xsize; % Col
mri.x_a = M(2,1)/mri.xsize;
mri.x_s = M(3,1)/mri.xsize;

mri.y_r = M(1,2)/mri.ysize; % Row
mri.y_a = M(2,2)/mri.ysize;
mri.y_s = M(3,2)/mri.ysize;

mri.z_r = M(1,3)/mri.zsize; % Slice
mri.z_a = M(2,3)/mri.zsize;
mri.z_s = M(3,3)/mri.zsize;

ic = [(mri.width)/2 (mri.height)/2 (mri.depth)/2 1]';
c = M*ic;
mri.c_r = c(1);
mri.c_a = c(2);
mri.c_s = c(3);
%--------------------------------------------------%

%-------- The stuff here is for convenience --------------

% 1-based vox2ras. Good for doing transforms in matlab
mri.vox2ras1 = vox2ras_0to1(M);

% Matrix of direction cosines
mri.Mdc = [M(1:3,1)/mri.xsize M(1:3,2)/mri.ysize M(1:3,3)/mri.zsize];

% Vector of voxel resolutions (Row-Col-Slice)
mri.volres = [mri.ysize mri.xsize mri.zsize];

% Have to swap rows and columns back
voldim = [mri.volsize(2) mri.volsize(1) mri.volsize(3)]; %[ncols nrows nslices]
volres = [mri.volres(2)  mri.volres(1)  mri.volres(3)];  %[dcol drow dslice]
mri.tkrvox2ras = vox2ras_tkreg(voldim,volres);


return;











