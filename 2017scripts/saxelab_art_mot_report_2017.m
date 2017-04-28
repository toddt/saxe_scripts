function artreport = saxelab_art_mot_report_2017(study,subjects,varargin)
% report = saxelab_art_mot_report_2017(study, subjects, bolds/model)
%
% This script reviews the art files of a subset of subjects in a study
% folder and compiles a report with all the relevant artifact and motion
% information.
%
% If using the default bolds option the script will look for wither all or
% a subset of the bold directories within each subjects' directory
% (depending on the bold input). If a model is specified, the relevant bold
% directories will be extracted from the SPM.mat filest
% The script then looks for the ART files in all the relevant bold 
% directories. If there are any missing art files or SPM files, the
% subject will be skipped and no art generation will be forced by this
% script. So no decision is going to be made and no new regressors will be
% constructed by this script. It only summarizes and reports what have been
% done in the process.
% 
% Note that art does not save the configurations used when it was run
% therefore this script can not provide any information on the
% configuration during the modeling part. So the user must know the
% thresholds and types provided to art at the time of modeling. Subjects
% that were analyzed using the new saxelab_generate_art_regressors script
% will have a config variable saved in the ART*.mat file so it can be
% checked later on per bold directory.
%
% All the artifacts data is compiled from the artifact regressors and the
% script has no knowledge for an artifact timepoint if it was originated
% from an intensity outlier or a motion outlier. 
%
% All motion data =is being computed from the motion regressors themselves
% and not from any output of arts. That is it is based on the data in the
% rp_rf0.txt file created by the original realignment process. 
%
% Inputs:
% study (char) - the study folder in either saxelab or saxelab2
% subjects (optional char or cell array) - allows to specify a subset of
%   the subjects in the study folder. if subjects is left out, the script
%   will search for all subjects starting with 'SAX_' in the study folder.
% bolds (numeric) - if the 3rd variable is left empty or does not exist, the
%   script will search for all the bold folders within the subject
%   directory and report on those (i.e. .../study/subject/bolds/0*).
%   bolds can be a numeric array, in which case the script will only look
%   for bold folders that match the numbers given. If there are multilple
%   subjects the script excepts either a cell array of bolds matching the
%   number of subjects or a single numeric array which will be assumed to
%   apply to all subjects.
%   If reporting on bold directories, no design relevant data will appear
%   in the report (artifacts by conditions).
% model (char) - if the 3rd variable exist and is a string, it will be 
%   treated as model.The script will search in each subject for a results 
%   folder named matching the specified model. (e.g. 'tom_results_normed') 
%
% Outputs:
% .csv file named 'art_mot_report_model_<date>.csv' will be created. If a 
% single subject was analyzed, it will be saved under subject/reports, if 
% multiple subjects, it will be saved in the study folder.
% report - a cell array that matches the report csv file will be returned.
% In the report, each row represent the data for a single subject's model.
% Each row has the following colomns of data:
% note: all translations are in mm and all rotation in degrees.
%
% subject - subject identifier
% model - the full name of the results folder 
% bold directories - the numbers of the bold directorues analyzed
% total artifacts - total number of artifact regressors in the model
% artifacts_by run - array of number of artifact regressors in each functional run
% mean translation - mean of the differential translation, meaned across runs and axis
% mean rotation - mean of the differential rotation, meaned across runs and axis
% mean distance - mean of the differential euclidean distance between time points
% mean angular rotation - mean of the differential 3d angle difference between timepoints using eauler's formula
% mean X - mean differential transition on single axis across runs
% mean Y - mean differential transition on single axis across runs 
% mean Z - mean differential transition on single axis across runs
% mean pitch - mean differential rotation on single axis across runs
% mean roll - mean differential rotation on single axis across runs
% mean yaw - mean differential rotation on single axis across runs
% mean trans by run - value per run, meaned across axis
% mean rot by run - value per run, meaned across axis
% mean distance by run - value per run
% mean angular rot by run - value per run
% artifacts by condition run # - for each condition in the design matrix,
%   gives the number of artifact timepoints dropped in that condition in 
%   the run. The last number is the number of rest timepoints dropped.
%   Condition is defined as the highest value in the design matrix for the
%   timepoint and rest is defined if all values in design matrix are under
%   0.1.



if nargin<1
    help saxelab_art_mot_report_2017;
    return;
end

if ~ischar(study)
    error('STUDY must be provided as a string input');
end

cwd = pwd;
%tmp = adir(sprintf('/mindhive/saxe*/%s',study));
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
% 

if ~isdir(study)
    error('Could not find study. Please provide full path.');
end


if nargin==1||isempty(subjects)
    % find all subjects which match.
%    search_string = fullfile(study,'SAX*');
%    subjects = adir(search_string);
    subjects = cellstr(spm_select('FPList',study,'dir','^SAX.*'));
%    if ~iscell(subjects)
    if isempty(subjects)
       error('Could not locate any subjects.');
    end
% elseif ischar(subjects) && any(strfind(subjects,'*'))
%     % a filter was given, find all subjects which match.
%     search_string = fullfile(study,subjects);
%     subjects = adir(search_string);
%     if ~iscell(subjects)
%         error('Could not locate any subjects.');
%     end
else
    if ischar(subjects)
        subjects = {subjects};
    end
    pot_subjects = subjects;
    subjects = {};
    if iscell(pot_subjects)
        for i = 1:length(pot_subjects)
%            tmp = adir(fullfile(study,pot_subjects{i}));
            tmp = cellstr(spm_select('FPList',study,'dir',['^' pot_subjects{i} '$']));
%            if ~iscell(tmp)
%                warning(sprintf('Could not find subject %s, will skip this subject\n',pot_subjects{i}));
            if isempty(tmp{1})
                error(['Could not find subject: ' pot_subjects{i}]);
            else
                subjects{end+1} = tmp{1};
            end
        end
        
    end
end
if length(subjects)<1
    error('Could not locate any matching subjects')
end

%This code is janky, but I'm not changing it now. --twt
modeled = 0;
if nargin > 2
    if ischar(varargin{1})&&~isempty(varargin{1})
        modeled = 1;
        model = varargin{1};
        for i = 1:length(subjects)
%            tmp = adir(fullfile(subjects{i},'results',model,'SPM.mat'));
            tmp = cellstr(spm_select('FPList',fullfile(subjects{i}, 'results', model), '^SPM.mat$'));
%            if iscell(tmp)
            if ~isempty(tmp{1})
                spms{i} = tmp{1};
                found_spm(i) = 1;
            else
%                 sprintf('could not find SPM.mat for subjects %s, will skip the subject\n',subjects{i});
%                 found_spm(i) = 0;
                error(['Could not find SPM.mat for subject: ' subjects{i}]);
            end
        end
%         subjects = subjects(logical(found_spm));
%         spms = spms(logical(found_spm));
    else
        bolds = varargin{1};
        if ~(isnumeric(bolds)||isempty(bolds)||iscell(bolds))
            error('the 3rd argument must be a model or bolds (numeric or empty for all bolds)')
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
    end
end
if exist('bolds','var')&&isempty(bolds{1})
    clear 'bolds'
end

fprintf('Found %i subjects.\n',length(subjects));
% use design in report
artreport = {'subject','model','bold dirs','total artifacts','artifacts_by run','mean translation', 'mean rotation',...
    'mean distance','mean angular rotation','mean X','mean Y','mean Z','mean pitch','mean roll','mean yaw',...
    'mean trans by run', 'mean rot by run','mean distance by run','mean angular rot by run'};
if modeled
    artreport = [artreport {'artifacts by condition run 1','artifacts by condition run 2','artifacts by condition run 3','artifacts by condition run 4'}];
end 

art_prefix = '';
for i = 1:length(subjects)
    close all 
    if modeled
        load(spms{i});
        tmp = regexp(spms{i},'/','split');
        study = fullfile('/',tmp{1:end-4});
        subject = tmp{end-3};
        model = tmp{end-1};
        artreport{i+1,1} = subject;
        artreport{i+1,2} = model;
        tmp = arrayfun(@(y) regexp(y,'/','split'),{SPM.xY.VY(:).fname},'UniformOutput',false);
        s_bolds = unique(cellfun(@(x) str2num(x{1}{end-1}),tmp));
        tmp = regexp(model,'with_','split');
        if length(tmp)>1
            art_prefix = regexp(tmp{2},'art','split');
            art_prefix = art_prefix{1};
        end
    else
        tmp = regexp(subjects{i},'/','split');
        artreport{i+1,1} = tmp{end};
        artreport{i+1,2} = 'not specified';
        model = 'no_model_specified';
        if ~exist('bolds','var')
%             tmp = adir(fullfile(subjects{i},'bold/0*'));
%             if isnumeric(tmp)
%                 fprintf('It seems that there are no bold folders in %s, skipping subject!',subjects{i});
%                 continue
%             end
            tmp = cellstr(spm_select('FPList',fullfile(subjects{i},'bold'),'dir','^0.*'));
            if isempty(tmp{1})
                error(['No bold folders for subject: ' subjects{i}]);
            end
            
            tmp = cellfun(@(x) regexp(x,'/','split'),tmp,'UniformOutput',false);
            s_bolds = unique(cellfun(@(x) str2num(x{end}),tmp));  
        else
            s_bolds = bolds{i};
        end
    end
    
    
    art_mot_files = {};
    for j = 1:length(s_bolds)
        %file = adir(fullfile(subjects{i},'bold',num2str(s_bolds(j),'%03.f'),[art_prefix 'art_regression_outliers_and_movement*.mat']));
        file = cellstr(spm_select('FPList', fullfile(subjects{i}, 'bold', num2str(s_bolds(j), '%03.f')), ...
        ['^' art_prefix 'art_regression_outliers_and_movement.*.mat']));
        %if iscell(file)
        if ~isempty(file{1})
            art_mot_files{end+1} = file{1};
        end
    end
    if length(art_mot_files)<length(s_bolds)
        fprintf('missing ART files for subject %s, skipping subject\n',subjects{i});
        artreport{i+1,3} = 'missing ART files';
        continue
    end

    arts_per_run = {};  % {conditions rest}
    mots_per_run = {};  % {x y z roll pitch yaw mean_trans mean_rotation mean_distance mean_auiler} 
    for j = 1:length(art_mot_files)
        mots = [];
        arts = [];
        arts_per_cond =[];
        load(art_mot_files{j})
        %first let's deal with artifacts
        arts = sum(R(:,1:end-6),2);
        if modeled
            conds = zeros(size(arts));
            cols = SPM.Sess(j).col(1:length(SPM.Sess(j).U));
            for tx = 1:length(SPM.Sess(j).row)
                if isempty(find(SPM.xX.X(SPM.Sess(j).row(tx),cols)>0.05))
                    conds(tx) = length(cols)+1;
                else
                    [tmp conds(tx)] = max(SPM.xX.X(SPM.Sess(j).row(tx),cols));
                end
            end
        else
            conds = ones(size(arts));
            cols = [];
        end
        arts = arts.*conds;
        for cx = 1:length(cols)+1
            arts_per_cond(cx) = length(find(arts==cx));
            images_per_cond(cx) = length(find(conds==cx));
        end
        images_per_run{j} = images_per_cond;
        arts_per_run{j} = arts_per_cond;
        % now let's deal with motion
        mots = diff(R(:,end-5:end),1,1);
        mots = [mots sqrt(sum(mots(:,1:3).^2,2))]; % adding differential euclidean distance
        mots = [mots acos((cos(mots(:,4)).*cos(mots(:,5)) + cos(mots(:,4)).*cos(mots(:,6)) + ...
        cos(mots(:,5)).*cos(mots(:,6)) + sin(mots(:,4)).*sin(mots(:,5)).*sin(mots(:,6)) - 1)/2)*180/pi]; % %use euler's formula to calculate a single rotation value 
        mots = nanmean(abs(mots),1);
        if isnan(mots)
            mots = NaN(1,8);
        end
        mots_per_run{j} = mots; 
    end
    artreport{i+1,3} = s_bolds;
    artreport{i+1,4} = sum(cellfun(@(x) sum(x), arts_per_run)); % total arts
    artreport{i+1,5} = cellfun(@(x) sum(x), arts_per_run); % total arts per run 
    artreport{i+1,6} = nanmean(cellfun(@(x) nanmean(x(1:3)), mots_per_run)); % mean translation across runs
    artreport{i+1,7} = nanmean(cellfun(@(x) mean(x(4:6))*180/pi, mots_per_run)); % mean rotation across runs
    artreport{i+1,8} = nanmean(cellfun(@(x) x(7), mots_per_run)); % mean distance across runs
    artreport{i+1,9} = nanmean(cellfun(@(x) x(8), mots_per_run)); % mean angular rotation across runs
    artreport{i+1,10} = nanmean(cellfun(@(x) x(1), mots_per_run)); % mean X translation
    artreport{i+1,11} = nanmean(cellfun(@(x) x(2), mots_per_run)); % mean Y translation
    artreport{i+1,12} = nanmean(cellfun(@(x) x(3), mots_per_run)); % mean Z translation
    artreport{i+1,13} = nanmean(cellfun(@(x) x(4)*180/pi, mots_per_run)); % mean roll
    artreport{i+1,14} = nanmean(cellfun(@(x) x(5)*180/pi, mots_per_run)); % mean pitch
    artreport{i+1,15} = nanmean(cellfun(@(x) x(6)*180/pi, mots_per_run)); %mean yaw
    artreport{i+1,16} = cellfun(@(x) nanmean(x(1:3)), mots_per_run); % mean trans per run
    artreport{i+1,17} = cellfun(@(x) nanmean(x(4:6))*180/pi, mots_per_run); % mean rotation per run
    artreport{i+1,18} = cellfun(@(x) x(7), mots_per_run); % mean distance per run
    artreport{i+1,19} = cellfun(@(x) x(8), mots_per_run); % mean angular per run
    if modeled
        for j = 1:length(arts_per_run)
            artreport{i+1,19+j} = arts_per_run{j}; % number of outlier per condition per run
        end
    end
end
if length(subjects) > 1
    mkdir(fullfile(study,'art_mot_reports'));
    f = fopen(fullfile(study,'art_mot_reports',['art_mot_report_' model '_' date '.csv']),'w');
else
    f = fopen(fullfile(subjects{1},'report',['art_mot_report_' model '_' date '.csv']),'w');
end
for i = 1:size(artreport,1)
    for j = 1:size(artreport,2)
        if ischar(artreport{i,j})
            fprintf(f,[artreport{i,j} ',']);
        else
            fprintf(f,'%G ',artreport{i,j});
            fprintf(f,',');
        end
    end
    fprintf(f,'\n');
end

    
fclose all;



