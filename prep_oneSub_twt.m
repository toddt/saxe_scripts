function study = prep_oneSub_twt(study)
% I have already converted dicoms for one subject. This script is to
% modify and implement the saxelab_prep_2014 script on the same subject.


%study.path = '/Users/hlk/Documents/MATLAB/saxelab_hlk/SUMM_hlk/SAX_summ_01/SAX_summ_01';
%study.subjects = {'SAX_new_01'};
study.prep = struct;
cd(study.path)

% I've decided it makes sense to use imbedded structures. That is, for the
% scripts that require their own variables that are not going to be used
% across scripts, it makes sense to have them in their own structure within
% the script strucutre. This is a way to have a simplified but more
% manageable "uconfig" variable.

% this script will also assume that the saxelab dicom script has been used
% to convert the dicom and therefore each subject should have the standard
% saxelab directories. If the dicoms were not converted using the saxelab
% dicom script, the user should make sure all of the standard dicom
% directories are present and/or run the makeSubjDirs.m script for each
% subject that she is using the prep script for.

% determine if the user has decided against using the default variables;
% the user will enter an alternative value for the following variables if
% she would like to use something other than the default. If she wants to
% use the default, she will not enter any of the other variables.

if ~isfield(study.prep,'realign')
    study.prep.realign = 1; % realign all images to first image of the run
end
if ~isfield(study.prep,'coreg')
    study.prep.coreg = 1; % coregister all images to first image of run
end
if ~isfield(study.prep,'normalize')
    study.prep.norm = 1; % normalize to MNI space
end
if ~isfield(study.prep,'smooth')
    study.prep.smooth = 1; % smooth using a fwhm mm smoothing kernel
end
if ~isfield(study.prep,'fwhm')
    study.prep.fwhm = 5; % the gaussian smoothing kernel to be used; value of fwhm in mm
end
fwhmFlag = 0;
fwhm = study.prep.fwhm;

% % if the data have already been pre-processed, need to process the dicoms
% % again.
% searchString = fullfile(subjects{1});
% tmp = dir(searchString);
% tmp = {tmp};
% if iscell(tmp) && ~isempty
%     dicom_oneSub_fx_hlk(study); % run the function to process the dicoms again but do not save the new study variable because you do not want it to interfere with the study variale already established in this script
% end

% now the preprocessing begins
subj = study.subjects{1};
cd(subj)
origDir = fullfile(study.path, subj, '3danat', 'orig'); % the directory that stores the the original anatomical images
repDir = fullfile(study.path, subj, 'report'); % the directory that stores reports generated by saxelab scripts
prepDate = strrep(datestr(clock),' ', '_'); % timestamp of preprocessing

% REALIGN
% find all of the bold images
bdirs = dir(['bold' filesep '0*']);
bdirs = {bdirs.name};
numBolds = length(bdirs);
bolds = cell(numBolds,1); % create a cell array to store bold image names
normBolds = bolds;
reslicBolds = bolds;
smoothNormBolds = bolds;
smoothReslicBolds = bolds;
for iBold = 1:numBolds
    boldDir = fullfile(study.path, study.subjects,'bold',bdirs{iBold});
    rawBolds = dir([boldDir{1} filesep 'f0*.img']);
    rawBolds = {rawBolds.name};
    bolds{iBold} = strcat(boldDir, '/',rawBolds,',1');
    if iBold == 1
        allFuncs = {bolds{iBold}{:}};
    else
        allFuncs = {allFuncs{:},bolds{iBold}{:}};
    end
end

% find the anatomical image
anatDir = fullfile(study.path,study.subjects,'3danat');
anat = dir([anatDir{1} filesep 's0*.img']); % find the name of the anatomical image
anat = {anat.name}; % create a string with the full path for the anatomical image
anat = strcat(anatDir,'/',anat,',1');

% start with estimating the realignment parameters for the functional images (realign the bolds to the first image)
if study.prep.realign
    estFlags = struct('quality',1,'rtm',0); % highest quality and don't register to mean
    for iEst = 1:numBolds
        %P = fullfile('bold',cellfun(@(x) x{end}, regexp(bolds{iEst},'/','split'),'UniformOutput',0));
        spm_realign(bolds{iEst},estFlags); % create the estimatation parameters for realignment --> file is saved as rp_*.txt
        spm_print(fullfile(repDir, sprintf('motionRun%s_%s', bdirs{iEst}(end-2:end), prepDate))); % save the motion report to the report directory
    end
end

% CO-REGISTER
% co-register the functional images to the anatomical images
% in saxelab_prep_2014.m if you are not normalizing, only the first image
% of the first run is registered to the anatomical; if you are normalizing,
% the first image of each run is registered to the anatomical. In this
% version, the first image of each run is registered to the anatomical
% regardless of if you are normalizing or not.
% so basically, what is happening here is the first image of each
% functional is registered to the anatomical then the remaining functional
% images are registered to the first functional of its respective run. The
% transformation parameters from the coregistration step are saved in the
% header information. These parameters can be viewed by running spm_vol and
% looking at the output from v.mat. The v.mat numbers change after the
% coregistration step but they are the same after the realignment step.
% coregistration.
if study.prep.coreg
    coreg.ref = {[anat{1} ',1']}; % the anatomical image
    for iCoreg = 1:numBolds
        coreg.source = {[bolds{iCoreg}{1} ',1']}; % this first functional image in the run
        coreg.other = bolds{iCoreg}(2:end); % the rest of the functional images in the run
        % saxelab_prep_2014.m uses spm_run_coreg_estimate but this function no
        % longer exists in spm12, it now uses spm_run_coreg.
        coreg.eoptions = spm_get_defaults('coreg.estimate'); % use the default options --> this variable needs to exist because if it's not, then spm_run_coreg will not work
        spm_run_coreg(coreg);
        spm_print(fullfile(repDir, sprintf('fun2anatCoreg_run%s_%s', num2str(bdirs{iCoreg}(end-2:end)), prepDate))); % save the realignment report for realigining the functionals to the anatmoical
    end
end

% NORMALIZE
% if the user is keeping the data in native space then we need to reslice
% and smooth. we might want to segment the native data too.
% if the user is transforming the data to a template, then we
% need to concatonate the coregistration transformation parameters and the
% normalization trasnformation parameters prior to slicing and smoothing.
% the normalization step in spm12 will concatonate coregistration and
% normalization information in addition to performing segmentation

% if the participant is in native space, need to realign, coregister,
% reslice, and smooth; if the participant is in nmi space, need to realign,
% coregister, normalize & slice, and smooth.

if ~study.prep.norm % if in standard space, slice
    for iReslice = 1:numBolds
        spm_reslice(bolds{iReslice});
    end
else % if normalized, normalize then concatonate coreg and norm parameters prior to slicing (includes segmentation)
    % first normalize the anatomical
    resample = {anat{1},allFuncs{:}};
    % set the variables for the job structure to feed to spm_run_norm
    norm.subj.vol            =   anat;      % anatomical path
    norm.subj.resample       =   resample'; % anatomical and functionals
    norm.eoptions.biasreg    =   0.0001;
    norm.eoptions.biasfwhm   =   60;
    norm.eoptions.tpm        =   {fullfile(spm('dir'),'tpm','TPM.nii')};
    norm.eoptions.affreg     =   'mni';
    norm.eoptions.reg        =   [0 1.0000e-03 0.5000 0.0500 0.2000];
    norm.eoptions.fwhm       =   0;
    norm.eoptions.samp       =   3;
    norm.woptions.bb         =   [-78 -112 -70; 78 76 85];
    norm.woptions.vox        =   [2 2 2];
    norm.woptions.interp     =   4;
    norm.woptions.prefix     =   'w';
    spm_run_norm_saxelab_2016(norm) % I'm going to try to run this under the assumption the spm will pull all the defaults.
end
% if we're going to delete the f0s for the sake of space, we should keep at
% one so we can look at it for QC purposes or we should take a picture of
% one

% SMOOTH
if study.prep.smooth
    if study.prep.norm
        addPref = 'w';
    else
        addPref = 'r';
    end
    smoothFuncs = char(allFuncs);
    for m = 1:length(allFuncs)
        curr = deblank(smoothFuncs(m,:));
        [path,name,ext] = fileparts(curr);
        inputImage = fullfile(path,[addPref name ext]);
        outputImage = fullfile(path, ['s' addPref name ext]);
        spm_smooth(inputImage,outputImage,fwhm);
    end
end

% GENMASK
% The next lines were copied from saxelab_prep_2014
sinfo = regexp(subj,'/','split');
system(sprintf('rm -rf %s',fullfile(subj,'mask')));
if study.prep.norm
    saxelab_genmask(study.path,subj);
else
    saxelab_genmask(study.path,subj,0);
end
prefix = 'f';
if ~ study.prep.norm prefix=['r' prefix]; end
if study.prep.norm prefix=['w' prefix]; end
if study.prep.smooth prefix=['s' prefix]; end
saxelab_checkreg(sinfo{end-1},sinfo{end},prefix);
spm_print(fullfile(repo_dir,sprintf('registration_check_%s',DOP)));

% DELETE UNNEEDED IMAGES
% also preserved from saxelab_prep_2014
if study.prep.realign
    system('rm -rf bold/*/f*');
end
if study.prep.norm||study.prep.smooth
    system('rm -rf bold/*/rf*');
end
if study.prep.norm&&study.prep.smooth
    system('rm -rf bold/*/wrf*');
end

end
