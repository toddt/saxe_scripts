function report = saxelab_generate_art_2017(study,subject,bolds,uconfig)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% report = saxelab_generate_art_2014(study,subjects,bolds,uconfig)
% 
% This script manages the creation of art regressors for subjects. It 
% should be run after preprocessing (one time or iteratively to decide on 
% desired thresholds and methog) and before modeling. The script will
% always result in 2 files saved in the bold directories the way that
% running art() usually does (i.e. art_regression_and_movement_outlier*.mat
% and art_regression_outliers*.mat) but with the difference that along with
% the R matrices there will be a config variable saved that has the
% configuration used when the regressors were generated. After generating
% and saving regressors, the script will call saxelab_art_mot_report and
% return the report compiled by that script for all the subjects that were
% analyzed.
%
% The user has a choice if to use art() as the way to pick outliers or if
% to use custom claculation of movement.
% Note that in the custom form, the calculations will be performed on the
% motion rp_* files and no global intensity outliers will be taken into
% account. The script will calculate euclidean distance between time points
% and 3d angular rotation and apply requested thresholds.
%
% All thresholds for either the art form or custom form are configurable.
%
% Inputs:
% study (char) - the study folder in either saxelab or saxelab2
% subjects (optional char or cell array) - allows to specify a subset of
%   the subjects in the study folder. if subjects is left out, the script
%   will search for all subjects starting with 'SAX_' in the study folder.
% bolds (optional numeric) - if left empty or does not exist, the
%   script will search for all the bold folders within the subject
%   directory and analyze those (i.e. .../study/subject/bolds/0*).
%   bolds can be a numeric array, in which case the script will only look
%   for bold folders that match the numbers given. If there are multilple
%   subjects the script excepts either a cell array of bolds matching the
%   number of subjects or a single numeric array which will be assumed to
%   apply to all subjects.
% uconfig - (optional cell array) - pairs of configuration and value thatn
% are given to the script to determine the value for those configurations
% i.e uconfig = {'diff_global',1,'global_threshold',2}
%
% What follows is a description of options that are currently configurable
% by the user (default value in parenthesis): 
% 
% % script general options:
%   overwrite (0)
%       > 1 or 0
%       > if 0 and standard art files already exist will prompt and ask if
%       to overwrite. if 1, will overwrite without asking
%   use_art (1)
%       > 1 or 0
%       > if 1, will use art() to pick outliers, if 0 will use the custom calculation
%   normed (1)
%       > 1 or 0
%       > if 1, will use normalized data, if 0, will use unnormalized data
%   custom_prefix ('')
%       > string
%       > custom prefix for searching for images (i.e., 'srf')
%   nestle_backward (0)
%       > number (TRs)
%       > number of TRs prior to identified outlier to regress out
%   nestle_forward (0)
%       > number (TRs)
%       > number of TRs after identified outlier to regress out
%   custom_name ('')
%       > string
%       > will add the string to the beginning of the name of all art files
%       (i.e. if custom_name is 'custom', art files will be named
%       'custom_art_regression...')
% art specific configurations (applicable if use_art = 1)
%   mask_type (1)
%       > 1 or 0
%       > if 1, will use genmask (see: saxelab_genmask)
%   global_mean (1)
%       > 1 or 0
%       > if 1, will use global mean as a criterion to find outliers
%   global_threshold (3)
%       > number (standard deviations)
%       > global signal threshold used to define outliers
%   use_diff_global (0)
%       > 1 or 0
%       > if 1, use global signal difference between scans to define
%         artifacts. if 0 looks at values as is in the motion files
%   composite_motion (1)
%       > 1 or 0
%       > if 1, use composite motion, if 0, look at single axis of
%         translation/rotation.
%   motion_threshold (2)
%       > number
%       > motion threshold used to define outliers. 
%         Note that in the case of composite motion this number will be the
%         threshold for composite motion (whatever that means) but if not
%         using composite motion this will be the translation threshold in
%         mm and the rotation threshold will be 0.02 radians which is hard
%         coded into art.
%   use_diff_motion (1)
%       > 1 or 0
%       > if 1, use motion difference between scans to define artifacts. if
%         0 looks at values as is in the motion files
% custom artifact specific configurations (applicable if use_art = 0)
%   translation_threshold (0.5)
%       > number (mm)
%       > differential distance threshold
%   rotation_threshold (0.5)
%       > number (degrees)
%       > differential 3d angular rotation threshold
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nir Jacoby 2/20/2014


if nargin<1
    help saxelab_generate_art_2017;
    return;
end

if iscell(subject)
    subject = char(subject);
end

cwd = pwd;
% tmp = adir(sprintf('/mindhive/saxe*/%s',study));
% if ~iscell(tmp)
%     error('Could not locate study.');
% elseif length(tmp) > 1
%     tmp = adir(sprintf('-t /mindhive/saxe*/%s',study));
%     study = tmp{1};
%     warning('Found more than one study with that name!')
%     warning(sprintf('Will analyze most recent one: %s',study));
% else
%     study = tmp{1};
% end

if ~isdir(study)
    error('Could not locate study; please provide full path.');
end


% if nargin==1||isempty(subjects)
%     % find all subjects which match.
%     search_string = fullfile(study,'SAX*');
%     subjects = adir(search_string);
%     if ~iscell(subjects)
%         error('Could not locate any subjects.');
%     end
%     fprintf('Found %i subjects to analyze.\n',length(subjects));
% elseif ischar(subjects) && any(strfind(subjects,'*'))
%     % a filter was given, find all subjects which match.
%     search_string = fullfile(study,subjects);
%     subjects = adir(search_string);
%     if ~iscell(subjects)
%         error('Could not locate any subjects.');
%     end
%     fprintf('Found %i subjects to analyze.\n',length(subjects));
% else
%     if ischar(subjects)
%         subjects = {subjects};
%     end
%     pot_subjects = subjects;
%     subjects = {};
%     if iscell(pot_subjects)
%         for i = 1:length(pot_subjects)
%             tmp = adir(fullfile(study,pot_subjects{i}));
%             if ~iscell(tmp)
%                 warning(sprintf('Could not find subject %s, will skip this subject\n',pot_subjects{i}));
%             else
%                 subjects{end+1} = tmp{1};
%             end
%         end
%         
%     end
% end
% if length(subjects)<1
%     error('Could not locate any matching subjects')
% end

if ~isdir(fullfile(study,subject))
    error('Could not locate subject.')
end


% Check that the bold directories if defined <bolds> are type DOUBLE or CELL
if ~exist('bolds','var')
    bolds = [];
end
% if ~isempty(bolds)
%     if ~isnumeric(bolds)&&~iscell(bolds)
%         error('BOLDS must be an array of DOUBLES or CELL array of arrays of DOUBLEs');
%     end
%     if isnumeric(bolds)
%         bolds = {bolds};
%     end
%     if iscell(bolds)
%         if length(bolds)==1
%             bolds = repmat(bolds,1,length(subjects));
%         elseif length(bolds)~=length(subjects)
%             error('Number of sets of bold directories specified does not match the number of subjects');
%         end
%     end
% end

if ~isempty(bolds)
    if ~isnumeric(bolds)
        error('If provided, BOLDS must be an array of doubles.');
    else
        bolds={bolds};
    end
end


%% Now set the configuration to default and update user configs.
config.overwrite = 0;
config.use_art = 1;
config.normed = 1;
config.custom_prefix = '';
config.nestle_backward = 0;
config.nestle_forward = 0;
config.custom_name = '';
% art configs
config.mask_type = 1;
config.global_mean = 1;
config.global_threshold = 3;
config.diff_global = 0;
config.composite_motion = 1;
config.motion_threshold = 2;
config.diff_motion = 1;
% custom outliers config
config.translation_threshold = 0.5;
config.rotation_threshold = 0.5;


if exist('uconfig','var')&&iscell(uconfig) % then they've included a custom configuration
    uco = fieldnames(config); % uco = user-configurable options
    for ovs = 1:2:length(uconfig)
        cov = uconfig([ovs,ovs+1]); % current configuration parameter-value pair
        dne = 0; % if 1, 'do not evaluate'
        if ~ischar(cov{1})
            warning('Extra variable %i is not a string',(ovs+1)/2);
            dne = 1;
        end
        % validation step 2: ensure that cov{1} is a variable that can be 
        % changed
        if ~dne&&~any(strcmp(cov{1},uco))
            warning('%s is not a user-configurable parameter',cov{1});
            dne = 1;
        end
        if ~dne
            config.(cov{1}) = cov{2};
        end
    end
end

fields = fieldnames(config);
if config.use_art
    for i = 15:16
        config = rmfield(config,fields{i});
    end
else
    for i = 8:14
        config = rmfield(config,fields{i});
    end
end
    

% if config.normed
%     prefix = 'swrf';
% else
%     prefix = 'srf';
% end
% if config.custom_prefix
%     prefix = config.custom_prefix;
% end

%FIX: Right now, everything is normed.
prefix = 'swrf';

%for s = 1:length(subjects)
    % look for bold folders
    if isempty(bolds)
%        bdirs = adir(fullfile(subjects{s},'bold','0*'));
        bdirs = cellstr(spm_select('FPList',fullfile(study,subject,'bold'),'dir','^0.*'));
%        if ~iscell(bdirs)
        if isempty(bdirs{1})
%             fprintf('could not locate any bold directories for %s, skipping subject.\n',subjects{s});
%             continue
            error('Could not locate any bold directories for subject.');
        end
    else
        bdirs = {};
%        for i = bolds{s}
        for i = bolds{1}
%            if iscell(adir(fullfile(subjects{s},'bold',sprintf('%03d',i))))
%                bdirs(end+1) = adir(fullfile(subjects{s},'bold',sprintf('%03d',i)));
%            end

            tmpDir = cellstr(spm_select('FPList',fullfile(study, subject, 'bold'), 'dir', sprintf('%03d',i)));
            if isempty(tmpDir{1})
                error('Could not find specified bold directories');
            else
                bdirs(end+1) = tmpDir{1};
            end
            
        end
%        if length(bdirs) < length(bolds{s})
%            fprintf('could not locate specified bold directories for %s, skipping subject.\n',subjects{s});
%            continue
%        end
    end
    
    % look for images
    images = {};
    for i = 1:length(bdirs)
%        tmp = adir(fullfile(bdirs{i},[prefix '*.img']));
        tmp = cellstr(spm_select('FPList',bdirs{i},['^' prefix '.*img']));
        
        %if iscell(tmp)
        if isempty(tmp{1})
            error('Could not locate images in all specified bold directories.')
        else
            images{end+1} = tmp;
        end
        
        %tmp = adir(fullfile(bdirs{i},['art_regression_outliers*' prefix '*.mat']));
        tmp = cellstr(spm_select('FPList',bdirs{i},['art_regression_outliers.*' prefix '.*.mat']));
        
        %if iscell(tmp)&~config.overwrite
        if ~isempty(tmp{1}) && ~ config.overwrite
%             tmp = questdlg('Old art files found! Delete?','art files found','Yes','No','No');
%             if ~strcmpi(tmp,'Yes');
%                 return;
%             end
            error('Old ART files found and overwrite not set. Stopping execution.');
        end
    end
%    if length(images) < length(bdirs)
%        fprintf('could not locate images in all specified bold directories for %s, skipping subject.\n',subjects{s});
%        continue
%    end

% If using art to generate art regressors
    if config.use_art
        % Check that the MASK is present
        if config.mask_type
            if config.normed
%                mask = adir(fullfile(subjects{s},'3danat','skull_strip_mask.img'));
                mask = cellstr(spm_select('FPList',fullfile(study,subject,'mask'),'^skull_strip_mask.nii'));
            else
                error('ART not currently supported for unnormed data. See TWT if you need this.');
                %mask = adir(fullfile(subjects{s},'3danat','unnormalized_skull_strip_mask.img'));
            end
%            if ~iscell(mask)
%                fprintf('could not find saxelab mask for %s, skipping subject.\n', subjects{s});
            if isempty(mask{1})
                error('Could not find saxelab mask for subject.');
            end
        else
            mask = -1;
        end
        
        % make an artifact directory
        
%        mkdir(fullfile(subjects{s},'art'));
        mkdir(fullfile(study,subject,'art'));
        if iscell(bolds)
%            task = num2str(bolds{s},'%02.f_');
            task = num2str(bold, '%02.f_');
            task = task(1:end-1);
        else
            task = 'all_bolds';
        end
%        mkdir(fullfile(subjects{s},'art',task));
%        cd(fullfile(subjects{s},'art',task));
        mkdir(fullfile(study, subject,'art',task));
        cd(fullfile(study, subject,'art',task));
        
        
        % construct a session file
        af = fopen('art_config001.cfg','w');
        fprintf(af,'sessions: %i\n',length(bdirs));
        fprintf(af,'global_mean: %i\n',config.global_mean);
        fprintf(af,'global_threshold: %.6f\n',config.global_threshold);
        fprintf(af,'motion_threshold: %.6f\n',config.motion_threshold);
        fprintf(af,'motion_file_type: 0\n');
        fprintf(af,'motion_fname_from_image_fname: 1\n');
        fprintf(af,'use_diff_motion: %i\n',config.diff_motion);
        fprintf(af,'use_diff_global: %i\n',config.diff_global);
        fprintf(af,'use_norms: %i\n',config.composite_motion);
        if iscell(mask)
            fprintf(af,'mask_file: %s\n',mask{1});
        end
%        fprintf(af,'output_dir: %s\n',fullfile(subjects{s},'art',task));
        fprintf(af,'output_dir: %s\n',fullfile(study, subject,'art',task));
        fprintf(af,'end\n');
        for i = 1:length(bdirs)
            fprintf(af,'session %i image ',i);
            fprintf(af,'%s ',images{i}{:});
            fprintf(af,'\n');
        end
        fprintf(af,'end\n');
        fclose(af);
        
        
        % run art 
%        fprintf('working on %s',subjects{s});
        fprintf('working on %s',subject);
        art('sess_file',fullfile(pwd,'art_config001.cfg'));
        close all
        for i = 1:length(bdirs)  % add config variable to saved art regressors .mat
%            artfiles = adir(fullfile(bdirs{i},'art_regression*.mat'));
            artfiles = cellstr(spm_select('FPList', bdirs{i}, '^art_regression.*.mat'));
            load(artfiles{1});
            nestle = 0;
            if config.nestle_backward||config.nestle_forward
                nestle = 1;
                mots = R(:,end-5:end);
                outliers = find(sum(R(:,1:end-6),2));
                for j = 1:config.nestle_backward
                outliers = union(outliers,outliers-1);
                end
                for j = 1:config.nestle_forward
                    outliers = union(outliers,outliers+1);
                end
                outliers = outliers(outliers>=1&outliers<=size(mots,1)); %make sure that nestling did not add values beyond the scope of timepoints
                R_nestled = zeros(size(mots,1),length(outliers));
                for j = 1:length(outliers)
                    R_nestled(outliers(j),j) = 1;
                end
            end
            if nestle
                R = [R_nestled mots];
            end           
            save(artfiles{1},'R','config');
            load(artfiles{2});
            if nestle
                R = R_nestled;
            end 
            save(artfiles{2},'R','config');
        end
    % If using custom motion outlier detection
    else 
        mot_files = {};
        for i = 1:length(bdirs)
%            tmp = adir(fullfile(bdirs{i},'rp_*.txt'));
%            if iscell(tmp)
            tmp = cellstr(spm_select('FPList',bdirs{i},'^rp_.*.txt'));
            if ~isempty(tmp{1})
                mot_files(end+1) = tmp(1);
            end
        end
        if length(mot_files) < length(bdirs)
%            fprintf('could not locate motion files in all specified bold directories for %s, skipping subject.\n',subjects{s});
%            continue
            error('Could not locate motion files in all specified bold directories.');
        end
        
        for i = 1:length(bdirs)
            mots = [];
            mots = load(mot_files{i});
            modiff = [zeros(1,6); diff(mots)]; %get between-tp motion
            translation = sqrt(sum(modiff(:,1:3).^2,2)); %compute distance
            rotation = acos((cos(modiff(:,4)).*cos(modiff(:,5)) + cos(modiff(:,4)).*cos(modiff(:,6)) + ...
            cos(modiff(:,5)).*cos(modiff(:,6)) + sin(modiff(:,4)).*sin(modiff(:,5)).*sin(modiff(:,6)) - 1)/2)*180/pi;
            outliers = unique([find(translation>config.translation_threshold); find(rotation>config.rotation_threshold)]); %compute 3d angle
            for j = 1:config.nestle_backward
                outliers = union(outliers,outliers-1);
            end
            for j = 1:config.nestle_forward
                outliers = union(outliers,outliers+1);
            end
            outliers = outliers(outliers>=1&outliers<=size(mots,1)); %make sure that nestling did not add values beyond the scope of timepoints
            R = zeros(size(mots,1),length(outliers));
            for j = 1:length(outliers)
                R(outliers(j),j) = 1;
            end
            tmp = regexp(images{i}{1},'/','split');
            tmp = regexp(tmp{end},'.img','split');
            tmp = tmp{1};
            save(fullfile(bdirs{i},['art_regression_outliers_' tmp]),'R','config');
            R = [R mots];
            save(fullfile(bdirs{i},['art_regression_outliers_and_movement_' tmp]),'R','config');
        end
    end
        
%end

cd(cwd);
%tmp = regexp(study,'/','split');
%study = tmp{end};
%subjects = cellfun(@(x) x{end}, regexp(subjects,'/','split'),'UniformOutput', 0);
try
%    report = saxelab_art_mot_report(study,subjects,bolds);
    artreport = saxelab_art_mot_report_2017(study,subject,bolds);
catch
    fprintf('All the art files were generated but reporting script failed. try running saxelab_art_mot_report with the same inputs\n')
end

if config.custom_name
%    rename_arts(study,subjects,bolds,config.custom_name,config.overwrite);
    rename_arts_2017(study,subject,bolds,config.custom_name,config.overwrite);

end

    
