function saxelab_model_2017(study,subjects,task,bolds,uconfig)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saxelab_model_2017(study,subjects,task,bolds,uconfig)
%
% This is the primary Saxelab batch function for creating and estimating a
% model. 
% 
% --------------------------BASIC FUNCTION--------------------------
% The script requires at least four arguments: 
% 
% saxelab_model_2017('study','subject','task',[bold])
%
% STUDY   - is a string, i.e., 'BLI'
% SUBJECT - is a string, i.e., 'SAX_BLI_01'
% TASK    - is a string, i.e., 'fba' such that there exist behavioral files
%           in the study's 'behavioural' directory, with the name 
%           <subject>.<task>.run#.mat
% BOLD    - is a numberic array, specifying the bold directories. 
%
% NOTE: similar to the old script, you can iterate over subjects by
% wrapping multiple subject strings in a CELL array (i.e., with curly
% braces). If this is the case, there must be either as many BOLD inputs as
% their are subject inputs (also wrapped in curly braces), or there must be
% exactly one--and the script will assume that each subject's analysis will
% use the same bold directories. 
%
% EXAMPLE ARGUMENTS:
%   saxelab_model_2017('BLI','SAX_BLI_01','fba',[5 9 14])
%   saxelab_model_2017('BLI',{'SAX_BLI_01','SAX_BLI_02','SAX_BLI_03'},'fba',{[5 9 14],[5 12 14],[7 9 11]})
%   saxelab_model_2017('BLI',{'SAX_BLI_01','SAX_BLI_02','SAX_BLI_03'},'fba',[5 9 14])
%   saxelab_model_2017('BLI',{'SAX_BLI_01','SAX_BLI_02','SAX_BLI_03'},'fba',[5 9 14],{'normed',0,'add_art',1})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------UCONFIG OPTIONS--------------------------
% The script supports a number of options which allow the analysis itself
% to be configured by the user. 
%
% These options are passed to the script as a series of parameter/value
% pairs in a cell array. It's not as hard as it sounds. The fifth argument
% to the function are where this occurs, and it consists of a series of
% such parameter/value pairs like:
%
%   {'param1',value,'param2',value,...,'paramN',value}
%
% For instance, if you're using unnormalized data and want to include
% artifact regressors, the fifth argument would look like:
%
%   {'normed',0,'add_art',1}
%
% and the entire function call would be
%
% saxelab_model_2017('BLI',{'SAX_BLI_01','SAX_BLI_02','SAX_BLI_03'},'fba',[5 9 14],{'normed',0,'add_art',1})
%
% What follows is a description of options that are currently configurable
% by the user (default value in parenthesis):
%
%   custom_name ('')
%       > string
%       > will name the results directory custom_name
%   clobber (0)
%       > 1 or 0
%       > if 1, will remove old results directory without prompting
%   rename_old_res (0)
%       > 1 or 0
%       > if 1, will rename old directories instead of overwriting them.
%   timing ('scans')
%       > string
%       > 'scans' or 'secs'
%       > switches between treating the timing information as TRs ('scans')
%         or seconds ('secs')
%   TR (2)
%       > number (seconds)
%       > the temporal resolution
%   content_delay (0)
%       > number (scans/seconds)
%       > introducing a delay to the expected content. This will shift the
%       onset and will be taken out of duration. If you know that your
%       mental content does not start until x TRs/sec into your stimuli, it
%       will make onset = onset+x and dur = dur-x because stimuli still
%       gets off screen at the anticipated time.
%   normed (1)
%       > 1 or 0
%       > if 1, will use normalized data, if 0, will use unnormalized data
%   custom_prefix ('')
%       > string
%       > custom prefix for searching for images (i.e., 'srf')
%   mask_type (1)
%       > 1 or 0
%       > if 1, will use genmask (see: saxelab_genmask)
%   filter_frequency (128)
%       > number (max seconds / cycle)
%       > highpass filter cutoff frequency
%   art_prefix ('')
%       > string
%       > custom prefix for art files (i.e., 'custom' will search for 
%       custom_art_regression...)
%   add_art (0)
%       > 1 or 0
%       > if 1, will include artifact regressors in the model
%   add_mot (0)
%       > 1 or 0
%       > if 1, will include motion regressors in the model
%   auto_corr_correct (0)
%       > 1 or 0
%       > if 1, will use an autoregressive model to remove unexplained
%         correlations
%   global_normalize (1)
%       > 1 or 0
%       > if 1, will account for global signal
%   inc_user_regs (0)
%       > 1 or 0
%       > if 1, will not attempt to include user regressors -- which can be
%         a problem in older studies where the user regressors field is not
%         properly specified. 
%   ignore_ips (0)
%       > 1 or 0
%       > if 1, will not validate the number of IPS on a run-by-run basis.
%         In some cases (like Bubbles3) the number of images in a run seems
%         to be at odds with the distribution of IPS for run, but the
%         images seems to sum to the correct number. Thus, by enabling
%         this, you can still model the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: Parse the initial inputs
% If there are no input arguments, display the help file
if nargin == 0
    help saxelab_model_2017
    return
end
% if there are less than 4 arguments, error out
if nargin<4
    disp('error: saxelab_model_2017 requires a minimum of 4 arguments');
    disp('''study'', ''subjects'', ''task'' & ''boldirs'' required');
    return
end
% Check that the study <study> is type CHAR
if iscell(study)
    if length(study)>1
        error('saxelab_model_2017 does not support iterating over studies');
    else
        study = study{1};
    end
elseif ~ischar(study)
    error('STUDY must be a STRING');
end
% Check that the subject <subjects> is type CHAR or CELL
if ~iscell(subjects)
    if ischar(subjects)
        subjects = {subjects};
    else
        error('subjects must be a STRING or CELL array of STRINGs');
    end
end
% Check that the task <task> is type CHAR
if iscell(task)
    if length(task)>1
        error('saxelab_model_2017 does not support iterating over tasks');
    else
        task = task{1};
    end
elseif ~ischar(task)
    error('TASK must be a STRING');
end
% Check that the bold directories <bolds> are type DOUBLE or CELL
if ~isnumeric(bolds)&&~iscell(bolds)
    error('BOLDS must be an array of DOUBLES or CELL array of arrays of DOUBLEs');
end
if isnumeric(bolds)
    bolds = {bolds};
end
if iscell(bolds)
    if length(bolds)==1
        bolds = repmat(bolds,1,length(subjects));
    elseif length(bolds)~=length(subjects)
        error('Number of sets of bold directories specified does not match the number of subjects');
    end
end

%% STEP 2.1: Serialize the processing of subjects
if length(subjects)>1
    for s = 1:length(subjects)
        if nargin>4
            saxelab_model_2017(study,subjects{s},task,bolds{s},uconfig);
        else
            saxelab_model_2017(study,subjects{s},task,bolds{s});
        end
    end
    return
end
tic;

%% STEP 2.2: Do some more variable processing
% Extract the variables
subj  = subjects{1};
bolds = bolds{1};
% make sure bolds is in the correct orientation, just in case
bolds = reshape(bolds,1,[]);

%% STEP 3: Initial check that everything is in order
% check that you can find the study
% study_dir = adir(['/mindhive/saxelab*/' study]);
%if ~iscell(study_dir)
if ~isdir(study)
    error('Could not locate study. Please provide full path.');
end
%study_dir = study_dir{1};

study_dir = study;

% check that you can find the subject
%subj_dir = adir(fullfile(study_dir,subj));
%if ~iscell(subj_dir)
subj_dir = fullfile(study_dir, subj);
if ~isdir(subj_dir)
    error('Could not locate subject');
end

% generate the results and report directory
mkdir(fullfile(subj_dir,'results'));
mkdir(fullfile(subj_dir,'report'));


% save the modeling command
comd = rlc();
save(fullfile(subj_dir,'report','model_command'),'comd');
% open the log


fclose all;
f=fopen(fullfile(subj_dir,'report',['saxelab_model_2017_' strrep(datestr(clock),' ','_') '.txt']),'w');
[~, user_name] = system('whoami');
fprintf(f,'Analysis by %s\n',user_name);
fprintf(f,comd); 

% check that you can find the bold directories
bdirs = {};
for i = bolds
%     if ~iscell(adir(fullfile(subj_dir,'bold',sprintf('%03d',i))))
%         error(sprintf('Could not locate bold directory %03d',i));
%     end
    tmp = cellstr(spm_select('FPList',fullfile(subj_dir,'bold'), 'dir', sprintf('^%03d$',i)));
    if isempty(tmp{1})
        error(sprintf('Could not locate bold directory %03d',i));
    else
%        bdirs(end+1) = adir(fullfile(subj_dir,'bold',sprintf('%03d',i)));
        bdirs(end + 1) = tmp;
    end
end


% check that you can find the same number of behavioral files as bold directories
% behavs = adir(fullfile(study_dir,'behavioural',[subj '.' task '.*.mat']));
% if ~iscell(behavs)
behavs = cellstr(spm_select('FPList', fullfile(study_dir,'behavioural'), ['^' subj '.' task '.*.mat']));
if isempty(behavs{1})
    error('Could not find any behavioral files');
elseif length(behavs)~=length(bolds)
    error('Number of bold directories does not match the number of behavioral files.');
end

% re-sort the behavioral to ensure they're in the correct order (an issue
% that can crop up when you have more than 9 runs.
[a__ b__] = sort(cellfun(@(x) str2num(x{end-1}),regexp(behavs,'\.','split')));
behavs = behavs(b__);

%% STEP 4: Parse remaining inputs
% define default values
def.custom_name         = '';   % permits you to specify a custom name (a string) as a results dir
def.clobber             = 0;    % if 1, will attempt to remove old results directories (if found)
def.rename_old_res      = 0;    % if 1, will attempt to rename old directories
def.timing              = 'scans'; % permits switching between treating the timing units as seconds or TRs ('scans')
def.TR                  = 2;    % TR in seconds
def.content_delay       = 0;    % how long after stimuli starts, relevant content starts
def.normed              = 1;    % Assume normalization
def.custom_prefix       = '';   % permits you to specify a custom prefix for the images that are being looked for (i.e., 'rf')
def.mask_type           = 1;    % if 1, use genmask. if 0, use SPM default
def.filter_frequency    = 128;  % default frequency (in seconds) for highpass filter
def.art_prefix          = '';   % will look for art files starting with prefix_art_regression...
def.add_art             = 0;    % if 1, will add artifact regressors
def.add_mot             = 0;    % if 1, will add motion regressors
def.add_cc              = 0;    % if 1, will add cc regressors
def.cc_num              = 5;    % number of principle components used as regs (specified in compcor)
def.auto_corr_correct   = 0;    % if 1, will correct for autocorrelations
def.global_normalize    = 1;    % if 1, will perform global intensity normalization by scaling
def.inc_user_regs       = 0;    % if 1, will include behavior-file specified user regressors
def.ignore_ips          = 0;    % if 1, will not perform by-run IPS validation

fprintf(f,'Study: %s\n',study);
fprintf(f,'Subject: %s\n',subj);
fprintf(f,'Task: %s\n',task);
fprintf(f,'Bolds: ');fprintf(f,'%03i ',bolds);fprintf(f,'\n');
fprintf(f,'Behavs: ');fprintf(f,'%s ',behavs{:});fprintf(f,'\n');
% user_config_vars is a list of user-configurable variables from the 
% command line--i.e., these are variables which are appropriate to modify
% by passing an executable string to saxelab_model_2017

% evaluate all strings passed to elaborate them into variable values
%   > if these do not work, throw a warning, but do not error out
user_config_vars    = fieldnames(def);
if nargin > 4
    for ovs = 1:2:length(uconfig)
        % uconfig is a series of param / value pairs
        cov = uconfig([ovs,ovs+1]); 
        dne = 0; % dne = 'do not evaluate'
        % validation step: ensure that cov{1} is a string
        if ~ischar(cov{1})
            warning('Extra variable %i is not a string',(ovs+1)/2);
            dne = 1;
        end
        % validation step 2: ensure that cov{1} is a variable that can be 
        % changed
        if ~dne&&~any(strcmp(cov{1},user_config_vars))
            warning('%s is not a user-configurable variable',cov{1});
            dne = 1;
        end
        if ~dne&&~isnumeric(cov{2})
            warning('%s is assigned a non-numeric value',cov{1});
        end
        if ~dne
            def.(cov{1}) = cov{2};
        end
    end
end

if ~any(strcmp({'scans','secs'},def.timing))
    warning('''%s'' is not an acceptable timing unit. Using ''scans''',def.timing);
    def.timing = 'scans';
end

% here we print more details about the configuration of the script to the
% report file...
for param = 1:length(user_config_vars)
    p_val = def.(user_config_vars{param}); % parameter value
    if isnumeric(p_val)
        if p_val==uint8(p_val)
            fprintf(f,'%s = %i\n',user_config_vars{param},p_val);
        else
            fprintf(f,'%s = %.4f\n',user_config_vars{param},p_val);
        end
    elseif ischar(p_val)
        if isempty(p_val)
            p_val = 'NONE/UNDEFINED';
            fprintf(f,'%s = %s\n',user_config_vars{param},p_val);
        else
            fprintf(f,'%s = %s\n',user_config_vars{param},p_val);
        end
    end
end


%% Step 5: Validate remaining requirements
% validate the behavioral file data against what's in the bold directories
images = {};
if def.normed
    prefix = 'swrf';
else
    prefix = 'srf';
end
if def.custom_prefix
    prefix = def.custom_prefix;
end
for b = 1:length(behavs)
    k = load(behavs{b});
    % if the user using the uconig otion to ignore ips, is this step
    % necessary?
    %if ~isfield(k,'ips')
     %   fprintf(f,'\n\nERROR: IPS not defined for behavioral %s',behavs{b});
      %  fclose(f);
       % error('IPS not defined for behavioral %s',behavs{b});
   % end
    if ~isfield(k,'spm_inputs')
        fprintf(f,'\n\nERROR: spm_inputs are not specified -- cannot proceed');
        fclose(f);
        error('spm_inputs are not specified -- cannot proceed');
    end
    comp_cons = 1;
    if ~isfield(k,'con_info')
        warning('con_info not specified -- will not attempt to generate contrasts');
        comp_cons = 0;
    end
%    images{b} = adir(fullfile(bdirs{b},[prefix '*.img']));
%    if ~iscell(images{b})
    images{b} = cellstr(spm_select('FPList',bdirs{b}, ['^' prefix '.*.img']));
    if isempty(images{b}{1})
        if def.normed
            fprintf(f,'\n\nERROR: data from run %s appear to be unnormalized',bdirs{b});
            fclose(f);
            error(sprintf('data from run %s appear to be unnormalized',bdirs{b}));
        else
            fprintf(f,'\n\nERROR: data from run %s appear to be missing',bdirs{b});
            fclose(f);
            error(sprintf('data from run %s appear to be missing',bdirs{b}));
        end
    end
    
    if length(images{b})~=k.ips&&~def.ignore_ips
        fprintf(f,'\n\nERROR: either too few or too many images in bold directory %s',bdirs{b});
        fclose(f);
        error(sprintf('either too few or too many images in bold directory %s',bdirs{b}));
    end
end

% Check that the MASK is present
if def.mask_type
    if def.normed
        mask = cellstr(spm_select('FPList',fullfile(study,subj,'mask'),'^skull_strip_mask.nii'));
        %mask = adir(fullfile(subj_dir,'3danat','skull_strip_mask.img'));
    else
        error('saxelab_model_2017 currently only supports normed data.');
        %mask = adir(fullfile(subj_dir,'3danat','unnormalized_skull_strip_mask.img'));
    end
    if ~iscell(mask)
        fprintf(f,'You want to use saxelab mask but it could not locate an appropriate mask. Run saxelab_genmask and model again.\n');
        error('You want to use saxelab mask but it could not locate an appropriate mask. Run saxelab_genmask and model again.\n');
    end
end


% Look for art and mot and cc files
art_mot_files = {};
art_files = {};
mot_files = {};
cc_files = {};


if isempty(def.art_prefix)
    art_prefix = '';
else
    art_prefix = [def.art_prefix '_'];
end
for i = 1:length(bdirs)
    if def.add_art
%        file_tmp = adir(fullfile(bdirs{i},[art_prefix 'art_regression_outliers_and_movement_' prefix '*.mat']));
%        if iscell(file_tmp)
        file_tmp = cellstr(spm_select('FPList',bdirs{i}, ['^' art_prefix 'art_regression_outliers_and_movement_' prefix '.*.mat']));
        if ~isempty(file_tmp{1})
            art_mot_files(i) = file_tmp(1);
        else
            fprintf(f,'\n\nERROR: you asked for artifact regressors to be modeled yet the matching art files seem to be missing');
            fclose(f);
            error('art files are missing. generate arts and try again.');
        end
        
%        file_tmp = adir(fullfile(bdirs{i},[art_prefix 'art_regression_outliers_' prefix '*.mat']));
%        if iscell(file_tmp)
        file_tmp = cellstr(spm_select('FPList',bdirs{i}, ['^' art_prefix 'art_regression_outliers_' prefix '.*.mat']));
        if ~isempty(file_tmp{1})
            art_files(i) = file_tmp(1);
        else
            fprintf(f,'\n\nERROR: you asked for artifact regressors to be modeled yet the matching art files seem to be missing');
            fclose(f);
            error('art files are missing. generate arts and try again.');
        end
    end
    
    
    if def.add_mot
%        file_tmp = adir(fullfile(bdirs{i},'rp_*.txt'));
%        if iscell(file_tmp)
        file_tmp = cellstr(spm_select('FPList',bdirs{i},'^rp_.*.txt'));
        if ~isempty(file_tmp{1})
            mot_files(i) = file_tmp(1);
        else
            fprintf(f,'\n\nERROR: you asked for motion regressors to be modeled but mot (rp_f*) files are missing');
            fclose(f);
            error('mot files are missing. generate arts and try again.');
        end
    end
    
    if def.add_cc
        if isempty(def.cc_num)
            cc_num = 5;
        else
            cc_num = def.cc_num;
        end
        
%         file_tmp = adir(fullfile(bdirs{i},[art_prefix num2str(cc_num) 'cc_*.txt']));
%         if iscell(file_tmp)
        file_tmp = cellstr(spm_select('FPList', bdirs{i}, ['^', art_prefix num2str(cc_num) 'cc_.*.txt']));
        if ~isempty(file_tmp{1})
            cc_files(i) = file_tmp(1);
        else
            fprintf(f,'\n\nERROR: you asked for cc regressors to be modeled but cc (cc_f*) files are missing');
            fclose(f);
            error('cc files are missing. generate arts and try again.');
        end
    end
    
end

%% Step 6: Specify a model directory
% construct a name for the model
model_dir = task;


if def.add_art && def.add_mot && def.add_cc
    model_dir = [model_dir '_with_' art_prefix 'art_and_mot_and_' num2str(cc_num) 'cc_reg'];
elseif def.add_art && def.add_cc
    model_dir = [model_dir '_with_' art_prefix 'art_and_' num2str(cc_num) 'cc_reg'];
elseif def.add_art && def.add_mot
    model_dir = [model_dir '_with_' art_prefix 'art_and_mot_reg'];
elseif def.add_mot && def.add_cc
    model_dir = [model_dir '_with_' num2str(cc_num) 'cc_and_mot_reg'];
elseif def.add_art
    model_dir = [model_dir '_with_' art_prefix 'art_reg'];
elseif def.add_cc
    model_dir = [model_dir '_with_' num2str(cc_num) 'cc_reg'];
elseif def.add_mot
    model_dir = [model_dir '_with_mot_reg'];
end


model_dir = [model_dir '_results'];
if def.normed
    model_dir = [model_dir '_normed'];
else
    model_dir = [model_dir '_unnormed'];
end
if ~isempty(def.custom_name)
    model_dir = def.custom_name;
end
% target dir
target_dir = fullfile(subj_dir,'results',model_dir);
if exist(target_dir,'dir')
    if def.rename_old_res
        % assign the current results directory the name
        % [target_dir '_old_' #] 
        % where # is the number of old results directories + 1
        sofar = adir([target_dir '_old_*']);
        if iscell(sofar)
            count = length(sofar)+1;
        else
            count = 1;
        end
        system(sprintf('mv %s %s',target_dir,[target_dir '_old_' sprintf('%02i',count)]));
        if exist(target_dir,'dir')
            fprintf(f,'\n\nERROR: Unable to relocate old results folder, bailing out to prevent overriding the old results directory!\n');
            fclose(f);
            fprintf('Unable to relocate old results folder, bailing out to prevent overriding the old results directory!\n');
            return
        end
    elseif def.clobber
        fprintf('Clobber is enabled, removing old directory...\n');
        system(sprintf('rm -rf %s',target_dir));
    else
        removeDir = questdlg('Old results directory found! Delete?','Dir with same name found','Yes','No','No');
        if strcmpi(removeDir,'Yes');
            system(sprintf('rm -rf %s',target_dir));
        else
            return;
        end
    end
end

% create the directory
mkdir(target_dir);

%% Step 7: Specify a model
fmri_spec.dir = {target_dir};
fmri_spec.timing.units = def.timing;
fmri_spec.timing.RT = def.TR;
fmri_spec.timing.fmri_t = 16; % number of time bins per scan
fmri_spec.timing.fmri_t0 = 8; % the bin in which the convolved hrf is sampled. 1 means sampling to the beginning of the image, 8 (t/2) is sampling to middle slice
fmri_spec.fact = struct('name', {}, 'levels', {});
fmri_spec.bases.hrf.derivs = [0 0];
fmri_spec.volt = 1; % OPTIONS: 1|2 = order of convolution; 1 = no Volterra
if def.global_normalize
    fmri_spec.global = 'Scaling';
else
    fmri_spec.global = 'None';
end
if def.auto_corr_correct
    fmri_spec.cvi = 'AR(1)';
else
    fmri_spec.cvi = 'none';
end
if def.mask_type
    spm_get_defaults('mask.thresh',-Inf); % set the mask threshold to -Inf to prevent it from thresholding
    fmri_spec.mask = {[mask{1}]};
end

% iterate over the behavioral files
for i = 1:length(behavs)
    k = load(behavs{i});
    fmri_spec.sess(i).scans = images{i};
    % iterate over the spm_inputs field
    for j = 1:length(k.spm_inputs)
        fmri_spec.sess(i).hpf = def.filter_frequency;
        fmri_spec.sess(i).cond(j).name = k.spm_inputs(j).name;
        fmri_spec.sess(i).cond(j).onset = k.spm_inputs(j).ons -1 + def.content_delay;   
        %the -1 fixes the fact that in saxelab convention onset of 1 meanf
        %first TR (stimuli starts at time 0) and in SPM, if the stimuli is
        %at the very first start onset should be 0.
        fmri_spec.sess(i).cond(j).duration = k.spm_inputs(j).dur - def.content_delay;
        fmri_spec.sess(i).cond(j).tmod = 0;
        fmri_spec.sess(i).cond(j).pmod = struct('name', {}, 'param', {}, 'poly', {});
        if isfield(k.spm_inputs(j),'pmod')&~isempty(k.spm_inputs(j).pmod)
            try
                fmri_spec.sess(i).cond(j).pmod(1).name = k.spm_inputs(j).pmod.name;
                fmri_spec.sess(i).cond(j).pmod(1).param = k.spm_inputs(j).pmod.param;
                try
                    fmri_spec.sess(i).cond(j).pmod(1).poly = k.spm_inputs(j).pmod.poly;
                catch
                    fmri_spec.sess(i).cond(j).pmod(1).poly = 1;
                end
            catch
                fmri_spec.sess(i).cond(j).pmod = struct('name', {}, 'param', {}, 'poly', {});
            end
        end
    end
    % iterate over the user_regressors field
    if isfield(k,'user_regressors')&&def.inc_user_regs
        for j = 1:length(k.user_regressors)
            try
                fmri_spec.sess(i).regress(j).name = k.user_regressors(j).name;
                fmri_spec.sess(i).regress(j).val = k.user_regressors(j).val;
            catch
                fmri_spec.sess = rmfield(fmri_spec.sess,'regress');
                break
            end
        end
    end
    % add motion / cc / artifacts regressors if needed
    if def.add_art
        k = load(art_files{i});
        artrs = k.R;
        if def.add_art
            if ~isfield(fmri_spec.sess(i),'regress')
                r_ind = 1;
            else
                r_ind = length(fmri_spec.sess(i).regress)+1;
            end
            for artr = 1:size(artrs,2)
                fmri_spec.sess(i).regress(r_ind).name = ['artifact_' num2str(artr) '_run_' num2str(i)];
                fmri_spec.sess(i).regress(r_ind).val = artrs(:,artr);
                r_ind = r_ind+1;
            end
        end
    end
    
    if def.add_cc
        ccrs = load(cc_files{i});
        
        cc_n = cell(1,cc_num);
        for numccs=1:cc_num
            %cc_n = {'cc1','cc2','cc3','cc4','cc5'};
            cc_n(1,numccs) = {['cc' num2str(numccs)]};
        end
        
        if ~isfield(fmri_spec.sess(i),'regress')
            r_ind = 1;
        else
            r_ind = length(fmri_spec.sess(i).regress)+1;
        end
        for ccr = 1:size(ccrs,2)
            fmri_spec.sess(i).regress(r_ind).name = [cc_n{ccr} '_run_' num2str(i)];
            fmri_spec.sess(i).regress(r_ind).val = ccrs(:,ccr);
            r_ind = r_ind+1;
        end
    end
    
    if def.add_mot
        motrs = load(mot_files{i});
        mot_n = {'x_trans','y_trans','z_trans','roll','pitch','yaw'};
        if ~isfield(fmri_spec.sess(i),'regress')
            r_ind = 1;
        else
            r_ind = length(fmri_spec.sess(i).regress)+1;
        end
        for motr = 1:size(motrs,2)
            fmri_spec.sess(i).regress(r_ind).name = [mot_n{motr} '_run_' num2str(i)];
            fmri_spec.sess(i).regress(r_ind).val = motrs(:,motr);
            r_ind = r_ind+1;
        end
    end
end

%% Step 8: Create the model specification
matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
save(fullfile(target_dir,'batch.mat'),'matlabbatch');
fprintf(f,'\nModel specification has begun. \n');
try
    spm_jobman('run',matlabbatch);
catch perror
    fprintf('\n\nERROR: Model specification has failed -- consult report.\n');
    fprintf(f,'\n\nERROR: Model specification has failed.\n');
    fprintf(f,'\n\nERROR Returned:\n%s\n',perror.message);
    fclose(f);
    return
end
cd(target_dir)

%% Step 9: Estimate the model
fprintf(f,'\nModel specification complete.\n');
load SPM.mat
fprintf(f,'\nModel estimation has begun. \n');
try
    spm_spm(SPM);
catch perror
    fprintf('\n\nERROR: Model estimation has failed -- consult report.\n');
    fprintf(f,'\n\nERROR: Model estimation has failed.\n');
    fprintf(f,'\n\nERROR Returned:\n%s\n',perror.message);
    fclose(f);
    return
end
fprintf('Model estimation complete.\n');

%% Step 10: Generate the contrasts
% in contrast to previous scripts, this script will not assume that 
% contrasts are the same between runs, instead it will be done on a run-by-
% run basis. 
if ~comp_cons
    fprintf('Not computing contrasts. Finished!\n');
    fprintf('Compute time: %02ih',floor(toc/3600));
    fprintf(' %02im',floor(rem(toc,3600)/60));
    fprintf(' %.2fs\n',rem(toc,60));
    fprintf(f,'\n\nNot computing contrasts (missing con_info). Finished!\n');
    fprintf(f,'Compute time: %02ih',floor(toc/3600));
    fprintf(f,' %02im',floor(rem(toc,3600)/60));
    fprintf(f,' %.2fs\n',rem(toc,60));
    fclose(f);
    return;
end

clear k
for b = 1:length(behavs)
    k(b) = load(behavs{b},'con_info'); 
    conlen = length(k(b).con_info);
    if b==1
        pconlen = conlen;
    end
    if pconlen~=conlen
        warning('There are an irregular number of contrasts across runs!');
        fprintf('Not computing contrasts. Finished!\n');
        fprintf('Compute time: %02ih',floor(toc/3600));
        fprintf(' %02im',floor(rem(toc,3600)/60));
        fprintf(' %.2fs\n',rem(toc,60));
        fprintf(f,'\n\nNot computing contrasts (irregular number of contrasts across runs). Finished!\n');
        fprintf(f,'Compute time: %02ih',floor(toc/3600));
        fprintf(f,' %02im',floor(rem(toc,3600)/60));
        fprintf(f,' %.2fs\n',rem(toc,60));
        fclose(f);
        return
    end
    pconlen= conlen;
end

load(fullfile(target_dir,'SPM.mat'));
% now that we've validated the contrast images, begin contrast estimation
for con = 1:length(k(1).con_info)
    cols = [];
    vals = [];
    c = zeros(1,size(SPM.xX.X,2));
    for run = 1:length(k)
        pos_ind = find(k(run).con_info(con).vals > 0);
        neg_ind = find(k(run).con_info(con).vals < 0);
        cols    = [cols SPM.Sess(run).col(pos_ind) SPM.Sess(run).col(neg_ind)];
        vals    = [vals k(run).con_info(con).vals(pos_ind) k(run).con_info(con).vals(neg_ind)];
    end
    c(cols) = vals;
    % make sure that the contrast is either balanced or main effect.
    csum = sum(c);
    if csum
        if any(c>0)&&any(c<0)
            c_orig = c;
            pval = find(c>0);cpv = abs(sum(c(pval)));
            nval = find(c<0);cnv = abs(sum(c(nval)));
            if csum > 1
                % scale down the positive values
                c(pval)=c(pval)./(sum(cpv)/sum(cnv));
            else
                c(nval)=c(nval)./(sum(cnv)/sum(cpv));
            end
            fprintf('Contrast is uneven! Adjusted from\n');
            fprintf('%5.2G ',c_orig);fprintf('Sum: %.2f',sum(c_orig));fprintf('\n');
            fprintf('to\n');
            fprintf('%5.2G ',c);fprintf('Sum: %.2f',sum(c));fprintf('\n');
        else
            fprintf('Contrast is uneven, but it looks like you''re\njust looking for a main effect of multiple conditions, so not gonna touch it.\n');
        end
    end
    if ~isfield(SPM,'xCon')||isempty(SPM.xCon)
        SPM.xCon = spm_FcUtil('Set', k(1).con_info(con).name, 'T', 'c', c',SPM.xX.xKXs);
    else
        SPM.xCon(end+1) = spm_FcUtil('Set', k(1).con_info(con).name, 'T', 'c', c',SPM.xX.xKXs);
    end
end
spm_contrasts(SPM);

fprintf('Computing contrasts complete.\n');
fprintf('Compute time: %02ih',floor(toc/3600));
fprintf(' %02im',floor(rem(toc,3600)/60));
fprintf(' %.2fs\n',rem(toc,60));
fclose all;
end