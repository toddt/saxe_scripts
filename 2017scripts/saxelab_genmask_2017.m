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
    threshold = 0.35;
end



% find the normalized anatomical image
anatDir = spm_select('FPList',pwd,'dir','^3danat.*');
anatImg = spm_select('FPList',anatDir,'^ws0.*img');
anatHdr = spm_select('FPList',anatDir,'^ws0.*hdr');
copyfile(anatImg, mskdir_name);
copyfile(anatHdr, mskdir_name);


% Find the warping information
deffield=spm_select('FPList',anatDir, '^y_s.*nii');    
copyfile(deffield, mskdir_name);
    
cd(mskdir_name);

%Use BET to skull strip normalized anatomical
fprintf('Creating skull-stripped normalized anatomical\n');
system(sprintf('/usr/local/fsl/bin/bet %s %s -R -f %.1f',anatImg,fullfile(pwd,'skull_strip_mask'),threshold));
system('gunzip *.nii.gz');

%Copy the preproc segmentations to the mask_dir as appropriate mask types
gray_matter = spm_select('FPList',anatDir, '^wc1s0.*nii');
white_matter = spm_select('FPList',anatDir, '^wc2s0.*nii');
CSF = spm_select('FPList',anatDir, '^wc3s0.*nii');
copyfile(gray_matter,fullfile(mskdir_name,'gray_matter_mask.nii'));
copyfile(white_matter,fullfile(mskdir_name,'white_matter_mask.nii'));
copyfile(CSF,fullfile(mskdir_name,'csf_mask.nii'));



