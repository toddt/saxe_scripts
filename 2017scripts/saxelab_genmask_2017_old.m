function saxelab_genmask_2017(study, subject, threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saxelab_getmask(study, subject[, threshold])
%
% generates mask using FSL's BET program. Also makes white matter, gray matter
% and CSF masks from the segmentation step of normalization
%
% Usage:
%
%   saxelab_genmask('/mindhive/saxelab2/TOM','SAX_TOM_01')
%   saxelab_genmask('/mindhive/saxelab2/TOM','SAX_TOM_01',0.4) (generates normalized mask
%   with a threshold of 0.4)
%

%% Find the study 
fprintf('Locating study & subject.\n');
root = spm_select('FPList',study,'dir',subject);

if isempty(root)
    error('Could not locate study/subject pair.');
elseif size(root,1) > 1
    error('Found more than one matching study/subject pair!')
end

mskdir_name = fullfile(root,'mask');
repo_dir = fullfile(root,'report');
mkdir(repo_dir);

cd(root);
fprintf('Removing old temporary masking directories...\n');
system(sprintf('rm -rfv %s',mskdir_name));
fprintf('Gathering files...\n');
mkdir(mskdir_name);

if ~exist('threshold','var')||isempty(threshold)    
    threshold = 0.5;
end

if normalized && ~iscell(adir('3danat/*seg*.mat'))
    error('Cannot produce normalized mask because normalization warp does not exist.');
end


% obtain the functional image of reference.
% If first run, create copies in orig folder, if not, take f0s from orig
% because we no longer keep f0s in study folders.

%bolds = adir('bold/*');bolds=bolds{1};
boldDirs = cellstr(spm_select('FPList','bold','dir','^0.*'));
bolds = boldDirs{1};

%if iscell(adir(fullfile(bolds,'f0*')))
funcs=cellstr(spm_select('FPList',bolds,'^f.*img'));
hdrs=cellstr(spm_select('FPList',bolds,'^f.*hdr'));
if ~isempty(funcs)
%    func = adir(fullfile(bolds,'f*.img'));func=func{1};
%    hdrs = adir(fullfile(bolds,'f*.hdr'));hdrs=hdrs{1};
    func1 = funcs{1};
    hdr1 = hdrs{1};
    system(['cp ' func1 ' 3danat/orig/']);
    system(['cp ' hdr1 ' 3danat/orig/']);
else 
    %FIX: This shouldn't happen right now, since we're not deleting f0s
    %from bold dirs or moving them to the 3danat orig dir, which was a
    %ridiculous place to put bolds in the first place.
    
%    func = adir('3danat/orig/f*.img');func=func{1};
%    hdrs = adir('3danat/orig/f*.hdr');hdrs=hdrs{1};
end

%copyfile(func,mskdir_name);
%copyfile(hdrs,mskdir_name);
copyfile(func1,mskdir_name);
copyfile(hdr1,mskdir_name);


% copy anatomical images of reference
%anat = adir('3danat/s0*.img');anat=anat{1};copyfile(anat,mskdir_name);
%hdrs = adir('3danat/s0*.hdr');hdrs=hdrs{1};copyfile(hdrs,mskdir_name);

% find the anatomical image
anatDir = spm_select('FPList',pwd,'dir','^3danat.*');
anatImg = spm_select('FPList',anatDir,'^s0.*img');
anatHdr = spm_select('FPList',anatDir,'^s0.*hdr');
copyfile(anatImg, mskdir_name);
copyfile(anatHdr, mskdir_name);

if normalized
    %    warp = adir('3danat/s*seg*mat');warp=warp{1};copyfile(warp,mskdir_name);
    
    % Find the warping information
    segmentation = spm_select('FPList',anatDir, '^s0.*seg8.mat');
    deffield=spm_select('FPList',anatDir, '^y_s.*nii');
    
    copyfile(segmentation,mskdir_name);
    copyfile(deffield, mskdir_name);
    
end

cd(mskdir_name);

%func = spm_vol(char(adir('f*.img')));

func=spm_select('FPList',mskdir_name,'^f.*img');
func=spm_vol(func);

fprintf('Creating skull-stripped anatomical\n');
%system(sprintf('/usr/local/fsl/bin/bet %s %s -R -f %.1f',char(adir('s0*.img')),fullfile(pwd,'ss[a]'),threshold));
system(sprintf('/usr/local/fsl/bin/bet %s %s -R -f %.1f',anatImg,fullfile(pwd,'ss[a]'),threshold));
system('gunzip *.nii.gz');

%anat = spm_vol(fullfile(pwd,'ss[a].nii'));
anat = spm_select('FPList',mskdir_name,'^ss\[a\].nii');
anat = spm_vol(anat);
%FIX: Stopped editing here. Maybe we want to copy the anatomical mask to
%the mask dir.

fprintf('Coregistering functional to anatomical.\n');
anon = @()spm_coreg(anat,func);
[a x] = evalc('anon()');
spm_print(fullfile(repo_dir,'coregistration_for_mask'));
anat = spm_vol(char(adir('s0*.img')));
M = inv(spm_matrix(x));
MM=spm_get_space(func.fname);
spm_get_space(func.fname,M*MM);

fprintf('Reslicing functional to anatomical.\n');
reslice_flags.which = 1;
reslice_flags.mean = 0;
spm_reslice({anat.fname;func.fname},reslice_flags);
realigned_func = spm_vol(char(adir('rf*.img')));
R = spm_read_vols(realigned_func);
realigned_func.fname = 'cor[f].img';
spm_write_vol(realigned_func,R);

fprintf('Binarizing coregistered functional.\n');
x = spm_vol(char(adir('cor[f].img')));
X = spm_read_vols(x);
X = X~=0;
x.fname = 'bin[cor[f]].img';
spm_write_vol(x,X);
%
fprintf('Validating labels...\n');
[tmp L] = bwlabeln(~X);
if L > 1
    warning('More than 1 zero-region present. You should check to mask to verify you do not have holes.');
end

fprintf('Binarizing skull-stripped anatomical...\n');
V = spm_vol(fullfile(pwd,'ss[a].nii'));
Vo = struct('fname',fullfile(pwd,'bin[ss[a]].img'),'dim',[V(1).dim(1:3)],'dt',[spm_type('float32'), 1],'mat',V(1).mat,'descrip','spm - algebra','mask',1);
Vo = spm_imcalc(V,Vo,'i1>0');

fprintf('Creating conjunction...\n');
Y = spm_read_vols(Vo);
new_mask = X&Y;
Vo.fname = 'unnormalized_skull_strip_mask.img';
Vo.descrip = 'Mask generated with saxelab_genmask';
cd ../3danat
spm_write_vol(Vo,new_mask);
cd ..
cd(mskdir_name);
if normalized
    Vo.fname = 'con[bin[ss[a]]_bin[cor[f]]].img';
    spm_write_vol(Vo,new_mask);
else
    cd(odir);
    return
end

fprintf('Applying warp to conjunction...\n');
p = load(char(adir('s*seg*mat')));
defaults=spm_get_defaults;
t = defaults.normalise.write;
t.prefix = 'n_';
t.interp = 0;
v0 = spm_write_sn('con[bin[ss[a]]_bin[cor[f]]].img',p,t);V = spm_read_vols(v0);
v0.fname = 'nrm[con[bin[ss[a]]_bin[cor[f]]]].img';spm_write_vol(v0,V);
cd ../3danat;
v0.fname = 'skull_strip_mask.img';spm_write_vol(v0,V);
t.interp = 1;
white_matter = adir(fullfile(root,'3danat','c2s*.img'));
grey_matter = adir(fullfile(root,'3danat','c1s*.img'));
fprintf('Applying warps to segmentation (grey matter)\n');
gmx = spm_write_sn(char(grey_matter),p,t);
fprintf('Binarizing grey matter into a mask...\n')
gmX = spm_read_Vols(gmx);
prob_thresh = 0.8;                      % probability requirement for inclusion in the masks
gmX = gmX >= prob_thresh; 
delete(gmx.fname);
gmx.fname = 'grey_matter_mask.img';
spm_write_vol(gmx,gmX);
fprintf('Applying warps to segmentation (white matter)\n');
gmx = spm_write_sn(char(white_matter),p,t);
fprintf('Binarizing white matter into a mask...\n')
gmX = spm_read_Vols(gmx);
gmX = gmX >= prob_thresh; 
delete(gmx.fname);
gmx.fname = 'white_matter_mask.img';
spm_write_vol(gmx,gmX);
t.interp = 0;
cd ..;
cd(mskdir_name);

fprintf('Applying additional warps for QC...\n');

v0 = spm_write_sn('ss[a].nii',p,t);V = spm_read_vols(v0);
v0.fname = 'nrm[ss[a]].nii';spm_write_vol(v0,V);

v0 = spm_write_sn('bin[cor[f]].img',p,t);V = spm_read_vols(v0);
v0.fname = 'nrm[bin[cor[f]]].img';spm_write_vol(v0,V);

v0 = spm_write_sn('cor[f].img',p,t);V = spm_read_vols(v0);
v0.fname = 'nrm[cor[f]].img';spm_write_vol(v0,V);

v0 = spm_write_sn(anat.fname,p,t);V = spm_read_vols(v0); 
v0.fname = 'nrm[a].img';spm_write_vol(v0,V);

system('rm -rf n_*');
fopen('complete_without_errors','w');
fclose all;
cd(odir);
end

