function saxelab_roi_picker_2016(study,subjects,model,con,roif,uconfig)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saxelab_roi_picker_2014(study,subjects,model,con,roif,uconfig)
%
%
% This is a batch script that picks ROIs from search spaces.
%
% The script will upload the approriate T-images and take following steps:
% 1. Threshold to match the requested p_val_thresh.
% 2. Mask the image with the search space
% 3. Eliminate clusters under cluster_size
% 4. Picks the highest T value in the remaining clusters as peak and the
%    cluster around it as the ROI
% 5. If required, mask the ROI with  a sphere_radius sphere around the peak 
%
% Inputs:
% all inputs can be left empty each empty variable will be prompted with a 
% selection GUI.
%
% study   - is a string, i.e., 'BLI'
% subjects - is a string/ cell array of subjects in the study folder or if
%          group_analysis = 1,folders from the RandomEffects folder
% model  - is a string - model name in the results folder, i.e.
%          'tomloc_with_art_reg_results_normed'
% con    - is a number - the number of contrast of interest 
% roif   - is a cell array with pointers to search spaces files (presumably
%          from the roi_library) 
% uconfig - cell array of option, value pairs {option#,value}
%
% Outputs:
% - subjects rois are saved under '<study>/<subject>/autoROI' with the name:
%       <ROI_name>_<model>_con<contrast number>_<suffix_>xyz.*
% - group rois are saved under '<study>/RandomEffects/autoROI/' with the name:
%      <ROI_name>_picked_<suffix_>xyz.* as both mat files and analyze image
% - subjects roi reports are saved under '<study>/<subject>/report/roi_picking' with the name:
%      ROI_defining_<date>.pdf
% - group roi reports are saved under '<study>/RandomEffects/report/roi_picking' with the name:
%      ROI_defining_<date>.pdf
% - log csv files will be saved in <study>/ROI with the name
%      saxelab_roi_picker_2014_log_<roi>_<sphere/cluster>_<date>.csv
%
% The script Constructs a job structure, then performs ROI picking. The job structure
% may be re-loaded into the function to repeat the picking.
% 
% job will consist of:
%   config - STRUCT, configuration parameters. 
%   study - CELL, containing study name
%   model - CELL, containing model name
%   rfx_folder - CELL, containing subject names / RFX directories
%   contrast - INT, the contrast number to pick from
%   conname - INT, the contrast name
%   rois - CELL, containing XYZmm, the ROI coordinates, for each ROI
%   roi_names - CELL of ROI names
%   roi_files - CELL of ROI files
%   roi_suffix - CHAR, full roi suffix (with _xyz, etc)
%   log_location - where to write the picker log to
%   log_fname - the log filename
%
% - the job file will be saved in <study>/autoROI_jobs/ with the name
%      roi_picking_job_<date>.m
% 
%
% options:
%
%   group_analysis
%       [boolean]
%       DEFAULT: 0
%       if 1 it permits picking from group ROIs. This will search for
%       "subjects" in the RandomEffects folder and omit the 'model' input.
%
%   p_val_thresh
%       [numeric]
%       DEFAULT: 0.001 (alpha)
%       p_val_thresh is the p-value threshold to be used in computing the
%       minimum T-value required for a voxel to be considered
%       suprathreshold. If you provide a number > 1, the script will assume
%       that it is the desired T threshold. 
%
%   clust_size
%       [numeric]
%       DEFAULT: 10 (voxels)
%       Specifies the minimum size of a contiguous cluster required for a
%       suprathreshold voxel to be included in the map.
%
%   op_order
%       [numeric]
%       DEFAULT: 1 
%       If 1,the order of things will be as described above. If 2 it will 
%       switch between steps 2 and 3. This will result in voxels
%       from outside the ROI being included in the peak.
%
%   old_style 
%       [boolean]
%       DEFAULT: 0
%       If 1, the script will pick the ROIs in the old style. If 0, then it
%       will pick ROIs in the current style.
%
%   whole_cluster
%       [boolean]
%       DEFAULT: 0
%       If 1, will select the whole cluster of the peak, as opposed to just
%       a sphere of the specified radius around it.
%
%   roi_radius
%       [numeric]
%       DEFAULT: 9 (mm)
%       Specifies the radius of the sphere to be drawn about the peak voxel
%       found in a given subject's contrast. 
%   
%   dilate_roi
%       [numeric]
%       DEFAULT: 1.00 
%       This must be a value between 1 and 0; if 1, it will take the ROI
%       mask as-is. If less than 1, will use the weighted ROI images (if it
%       can find them) to dilate the ROI, taking only the top
%       (dilate_roi*100)% of voxels. 
%
%   custom_suffix
%       [char]
%       DEFAULT: ''
%       Defining a custom suffix allows the user to add a custom string to
%       all ROIs picked in the current session. This is useful if there
%       will be multiple iterations of ROI picking, perhaps from the same
%       model but using different contrasts. Or the user may add in a date
%       of ROI picking, for instance. The custom suffix is inserted before
%       the _xyz.mat but after all other parts of the filename, or at the
%       very end if the _xyz suffix is omitted (see incl_xyz_stuff)
%
%   render 
%       [boolean]
%       DEFAULT: 1
%       If 1, will render images of the ROIs as they are produced, as well
%       as saving them in <subject/RFX_directory>/report/*, as PDF files.
%
%   top_voxels
%       [numeric]
%       DEFAULT: 0
%       If value is X(>0) the script will pick the top X voxels in the
%       search space. T will be thresholded at 0 and voxels don't have tuo
%       be contiguous
%-------------------------------------------------------------------------%
% Adapted from autoROI_picker
% Nir Jacoby 2/21/2014

cwd             = pwd;
close all force;

%% Check if the function is to run an existing job or create a new one.
% Case input is just a .mat file containing an saxelab_roi_picker_2014 job struct
if nargin == 1 && ischar(study) && exist(study,'file') && ~exist(study,'dir')
    fprintf('Loading job\n');
    load(study);
% Case input is just an saxelab_roi_picker_2014 job struct
elseif nargin == 1 && isstruct(study)
    fprintf('Job accepted\n');
    job = study;
% Case the input is the string 'load', opens a SPM selection Dialogue to
% select the .mat file containing the job.
elseif nargin == 1 && strcmpi(study,'Load')
    load(spm_select(1,'mat','Select job'));
else
%% Constract a new job
%% Config setting
% First set all configs to default value.
    config.group_analysis   = 0;
    config.p_val_thresh     = 0.001;
    config.clust_size       = 10;
    config.op_order         = 1;
    config.old_style        = 0;
    config.whole_cluster    = 0;
    config.roi_radius       = 9;
    config.dilate_roi       = 1;
    config.custom_suffix    = '';
    config.render           = 1;
    config.top_voxels       = 0;
    
% Check if there's a uconfig cell array (for costume config pairs)
    if exist('uconfig','var')&&~iscell(uconfig)
        fprintf('\n\tWARNING: Your user configuration is not a cell array!\n');
    end
% Case there is uconfig, check that it obeys the config pairs structure and
% all the fields refer to available configs.
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
    
    if config.top_voxels
        config.p_val_thresh = 0.5;    %	T threshold at 0
        config.whole_cluster = 1;         %no sphere masking
        config.clust_size = 1;      %no size thresholding
    end
        
% Write the configs into the job struct
    job.config = config;
    
%% Choose Study
    if (~nargin||(exist('study','var')&&isempty(study)))
        % no study has been specified
        study = inputdlg('What is the study name?');
        if isempty(study)
            study = '*';
        end

        pos_stu = adir(fullfile('/mindhive','saxelab*',['*' study{:} '*']));
        if ~iscell(pos_stu)
            pos_stu = adir(fullfile('/*','saxelab*',study{:}));
            if ~iscell(pos_stu)
                error('Could not locate study.\n');
            end
        end
        if length(pos_stu) > 1
            % the user must disambiguate which study they wish
            study = pos_stu(listdlg('ListString',pos_stu,'SelectionMode','single','PromptString','Multiple studies found, please choose one','ListSize',[500,300]));
            if isempty(study)
                error('No study chosen.\n');
            end
        else
            study = pos_stu;
        end
    else
        if ~iscell(study)
            study = {study};
        end

        pos_stu = adir(fullfile('/mindhive','saxelab*',['*' study{:} '*']));
        if ~iscell(pos_stu)
            pos_stu = adir(fullfile('/*','saxelab*',study{:}));
            if ~iscell(pos_stu)
                error('Could not locate study.\n');
            end
        end
        if length(pos_stu) > 1
            % the user must disambiguate which study they wish
            study = pos_stu(listdlg('ListString',pos_stu,'SelectionMode','single','PromptString','Multiple studies found, please choose one','ListSize',[500,300]));
            if isempty(study)
                error('No study chosen.\n');
            end
        else
            study = pos_stu;
        end
    end
    job.study = study;
    
% make sure that all non existing input variables exist and empty
    cd(study{1});
    if ~exist('subjects','var')   
        subjects = '';
    end
    if ~exist('model','var')
        model = '';
    end
    if ~exist('con','var')
        con = [];
    end
    if ~exist('roif','var')
        roif = {};
    end
    
 %% Choose subjects and model/ Group Analysis    
    if config.group_analysis
        model = {''};
        pos_sub = adir(['RandomEffects/' subjects '*']);
        if ~iscell(pos_sub)
            error('Could not locate RFX folders\n');
        elseif length(pos_sub) > 1
            subjects = pos_sub(listdlg('ListString',pos_sub,'SelectionMode','multiple','PromptString','Multiple RFX folders found, please select desired folders','ListSize',[300,300]));
        else
            subjects = pos_sub;
        end
    else
        if ~iscell(subjects)
            pos_sub = adir([subjects '*/results/*/SPM.mat']);
            pos_sub = unique(cellfun(@(x) x{1}, regexp(pos_sub,'/','split'),'UniformOutput',false));
            if ~iscell(pos_sub)
                error('Could not locate subject(s)\n');
            elseif length(pos_sub) > 1
                subjects = pos_sub(listdlg('ListString',pos_sub,'SelectionMode','multiple','PromptString','Multiple subjects found, please choose desired subjects','ListSize',[300,300]));
            else
                subjects = pos_sub;
            end
        end
 %       model = strrep(model,'/','');
        pos_model = adir([subjects{1} '/results/*' model '*']);
        pos_model = unique(cellfun(@(x) x{end},regexp(pos_model,'/','split'),'UniformOutput',0));
        if ~iscell(pos_model)
            error('Could not locate model\n');
        elseif length(pos_model) > 1
            model = pos_model(listdlg('ListString',pos_model,'SelectionMode','single','PromptString','Please select the model','ListSize',[500,300]));
        else
            model = pos_model;
        end
% check which of the subjects were modeled for the desired model and choose
% only those as the subjects for the job
        missing = [];
        notmiss = [];
        for i = 1:length(subjects)
            if ~exist(fullfile(subjects{i},'results',model{1},'SPM.mat'))
                missing(end+1) = i;
            else
                notmiss(end+1) = i;
            end
        end
        if ~isempty(missing)
            fprintf('WARNING: Could not find results directories for the following subjects:\n');
            for i = missing
                fprintf('\t%s\n',subjects{i});
            end
            fprintf('They will be excluded.\n');
            subjects = subjects(notmiss);
        end
    end
    job.model = model;
    job.subjects = subjects;
%% Choose contrast    
%load the SPM.mat for the desired job (either the RFX or a subject SPM.mat) to see available contrasts 
    if ~config.group_analysis
        load(fullfile(subjects{1},'results',model{1},'SPM.mat'));
    else
        load(fullfile(subjects{1},model{1},'SPM.mat'));
    end
    if isempty(con)
        contrast = listdlg('ListString',{SPM.xCon.name},'SelectionMode','single','PromptString','Please select the contrast','ListSize',[500,300]);
    elseif ~isnumeric(con)
        error('Contrast number must be an integer!\n');
    else
        contrast = con;
    end
    conname = SPM.xCon(contrast).name;
    job.contrast = contrast;
    job.conname = conname;
    
%% Choose hypotheses space ROIs
    if isempty(roif)
        roif = cellstr(spm_select([1 Inf],'mat','Select ROI xyz files','','/mindhive/saxelab/roi_library/'));
    end
    if isempty(roif{1})
        error('No ROIs selected!');
    end
    rois = {};
    names = {};
    for i = 1:length(roif)
        tmp = load(roif{i});
        tmpf = fieldnames(tmp);
% some hypotheses ROI files contain the xY struct and others don't, depends
% on the ROI library/ file that is used
        if any(strcmpi('xY',tmpf))  
            tmp4 = load(roif{i},'xY');
            rois{end+1} = tmp4.xY;
            [null name null] = fileparts(roif{i});
            name = strrep(name,'_xyz','');
            names{end+1} = name;
        elseif any(strcmpi('roi_xyz',tmpf))
            tmp4 = load(roif{i},'roi_xyz');
            if size(tmp4.roi_xyz,1)>3
                tmp4.xY.XYZmm = tmp4.roi_xyz';
            else
                tmp4.xY.XYZmm = tmp4.roi_xyz;
            end
            rois{end+1} = tmp4.xY;
            [null name null] = fileparts(roif{i});
            name = strrep(name,'_xyz','');
            names{end+1} = name;
        else
            fprintf('Your hypothesis spaces maybe be ''old-style'', and are not able to be used\n');
            error('Unrecognized ROI type!');
        end
    end
    job.rois = rois;
    job.roi_names = names;
    job.roi_files = roif;
% determine if there exists weighted rois, which are the same as ROIs,
% only with a 'wxyz' intead of 'xyz'
    if config.dilate_roi < 1
        wfiles = {};
        for rfile = 1:length(job.roi_files)
            null = regexp(job.roi_files{rfile},'_xyz.mat','split');
            wfile{rfile} = [null{1} '_wxyz.img'];
            if ~exist(wfile{rfile},'file')
                fprintf('\n\tWarning: Could not find weighted file for some of the ROIs, skipping dilation if requested\n');
                wfile = {};
                break
            end
        end
        if ~isempty(wfile)
            fprintf('Dilating ROIs\n');
            for rfile = 1:length(job.roi_files)
                wr = spm_vol(wfile{rfile});
                [WR XYZ] = spm_read_vols(wr);
                tthresh = prctile(WR(~~WR),100-(config.dilate_roi*100)); % tthresh is the minimum t-value threshold to dilate the ROI appropriately
                rois{rfile} = XYZ(:,find(WR)>=tthresh);
            end
            job.rois = rois;
            job.roi_wfiles = wfile;
        else
            fprintf('While ROI dilation was requested, the weighted ROI files could not be located and so this step is being skipped.\n');
        end
    end
%% Assign location and suffix for the resulting ROI files and create the directory where the results will be created
    if ~config.group_analysis
        roi_suffix = ['_' job.model{1} '_' sprintf('con%02i',job.contrast) '_'];
        if config.custom_suffix
            roi_suffix = [roi_suffix config.custom_suffix '_'];
        end
        roi_suffix = [roi_suffix 'xyz'];
    else
        if config.custom_suffix
            roi_suffix = ['_' config.custom_suffix '_xyz'];
        else
            roi_suffix = '_picked_xyz';
        end
    end
    job.roi_suffix = roi_suffix;
    job.log_location = fullfile(job.study{1},'ROI');
    if ~exist(job.log_location,'dir')
        mkdir(job.log_location);
    end
    if ~isempty(job.roi_names)
        for i = 1:length(job.roi_names)
            job.log_fname{i} = ['saxelab_roi_picker_2016_log_' job.roi_names{i} '_' job.model{1}];
            if config.top_voxels
                job.log_fname{i} = [job.log_fname{i} '_top_voxels' date];
            elseif config.whole_cluster
                job.log_fname{i} = [job.log_fname{i} '_cluster' date];
            else
                job.log_fname{i} = [job.log_fname{i} '_sphere' date];
            end
            if ~isempty(job.config.custom_suffix)
                job.log_fname{i} = [job.log_fname{i} '_' job.config.custom_suffix];
            end
            job.log_fname{i} = [job.log_fname{i} '.csv'];
            if exist(fullfile(job.log_location,job.log_fname{i}),'file'),
                f=fopen(fullfile(job.log_location,job.log_fname{i}),'a');
            else
                f=fopen(fullfile(job.log_location,job.log_fname{i}),'w');
            end
            if ~config.whole_cluster
                fprintf(f,'ROIs chosen for %s contrast %s at p = %G with %Gmm radius sphere\n',job.roi_names{i},job.conname,job.config.p_val_thresh,job.config.roi_radius);
            elseif config.top_voxels
                fprintf(f,'ROIs chosen for %s contrast %s at p = %G picking top %i voxels\n',job.roi_names{i},job.conname,job.config.p_val_thresh,job.config.top_voxels);
            else
                fprintf(f,'ROIs chosen for %s contrast %s at p = %G from a whole cluster\n',job.roi_names{i},job.conname,job.config.p_val_thresh);
            end
            if config.top_voxels
                fprintf(f,'subject,peak x,peak y,peak z,n voxels,mean T value,max T value,min Tvalue\n');
            else
                fprintf(f,'subject,peak x,peak y,peak z,n voxels,max T value\n');
            end
            fclose(f);
        end
    end
    mkdir(fullfile(study{1},'autoROI_jobs'));
    save(fullfile(study{1},'autoROI_jobs',['roi_picking_job_' strrep(strrep(strrep(datestr(clock),' ','_'),':','_'),'-','_') '.mat']),'job');
end

%% Process the job - Pick ROIs for each subject/ the group
cur_count = 0;
for s = 1:length(job.subjects)
    job = pick_ROIs(job.study{1},job.subjects{s},job.model{1}, ...
                       job.contrast,job.rois,job.roi_names,job.config, ...
                       job.roi_suffix,cur_count,job);
    cur_count = cur_count + length(job.rois);
end
cd(cwd);
end

%% The ROI picking procedure
function job = pick_ROIs(study,subject,model,contrast,rois,names,config,roi_suffix,cur_count,job)
%% Prepeartion activities before processing
% create the directory for the ROIs and navigate to the folder where the model results are. 

cd(study);
cd(subject)
mkdir('autoROI');
targ = fullfile(pwd,'autoROI');
if ~config.group_analysis
    cd(fullfile('results',model));
end
% report to command window which subject is being processed
tmp = regexp(subject,'/','split');
if ~config.group_analysis
    try
        subject = tmp{find(cellfun(@(x) ~isempty(x),tmp),1,'last')};
    catch
        subject = tmp{1}; 
    end
else
    subject = strrep(subject,[study '/'],'');
    subject = strrep(subject,study,'');
end
fprintf('%s...\n',subject);
% Create folder for reports and name for report file
pfilename = fullfile(study,subject,'report');mkdir(pfilename);
pfilename = fullfile(pfilename,'ROI_picking');mkdir(pfilename);
if config.top_voxels
    filename = fullfile(pfilename,sprintf(['ROI_top%Gvoxels_picking_' date],config.top_voxels));
elseif config.whole_cluster
    filename = fullfile(pfilename,['ROI_cluster_picking' date]);
else
    filename = fullfile(pfilename,sprintf(['ROI_%Gmm_sphere_picking_' date],config.roi_radius));
end
%% Loading images and prepare ROIs
% load up the T-image
if ~config.group_analysis
    t = spm_vol(fullfile(study,subject,'results',model,sprintf('spmT_%04i.img',contrast)));
else
    t = spm_vol(fullfile(study,subject,sprintf('spmT_%04i.img',contrast)));
end
[T XYZ] = spm_read_vols(t);

% threshold the spmT image
if config.p_val_thresh < 1  %is p val > 1 assume you were provided with T value.
    null = regexp(t.descrip,'[[]]','split');
    t_crit = tinv(1-config.p_val_thresh,str2num(null{2})); %set the T value threshold corresponding to the p thresh chosen
else
    t_crit = config.p_val_thresh;
end
if ~config.top_voxels
    T(~isfinite(T))=0; %Turns all NaNs to 0
    T(T<t_crit)=0; %Turns all values under the threshold to 0
else
    T(~isfinite(T))=-999;
end


if isempty(find(T))  %check that there were any voxels above T-Threshold
    fprintf('No suprathreshold voxels!\n');
    return
end

% transform the ROIs from millimeter space into index space
MAT    = t.mat;
IMAT   = inv(MAT);
for r = 1:length(rois)
    null = rois{r};null=null.XYZmm;
    null(4,:) = 1;
    null = IMAT*null;
    null = sub2ind(t.dim,null(1,:),null(2,:),null(3,:));
    roi_inds{r} = null;
end

%% Process peaks and ROI picking according to the requested order of actions
if config.op_order == 2||config.old_style
    % eliminate clusters smaller than the critical cluster size NOW, prior
    % to masking. 
    L = bwlabeln(~~T);
    clust_IDs = unique(L(~~L));
    for cID = clust_IDs'
        if numel(find(L==cID))<config.clust_size
            L(L==cID)=0;
        end
    end
    T = T.*~~L;
    if isempty(find(T))
        fprintf('No sufficiently sized clusters!\n'); %in all the image, regardless of ROIs
        return
    end
end

% find all the peaks, in case you need them later
peaks_ind = find(imregionalmax(T));
[px py pz] = ind2sub(size(T),peaks_ind);
peaks_xyz = [px py pz];

for r = 1:length(rois)
    % begin iterating through the rois
    stop = 0;
    cur_count = cur_count+1;
    fprintf('\tROI: %s\n',names{r});
    % create an ROI mask
    RM = zeros(t.dim);
    RM(roi_inds{r}) = 1;
    mT = T.*RM;
    if isempty(find(mT))
        fprintf('No suprathreshold voxels within ROI!\n');
        stop = 1;
    end
    % mT is now the ROI masked T-image
    

    % Check that the cluster size is larger than min cluster size required (if regular order) 
    if config.op_order == 1&&~config.old_style
        L = bwlabeln(~~mT);
        clust_IDs = unique(L(~~L));
        for cID = clust_IDs'
            if numel(find(L==cID))<config.clust_size
                L(L==cID)=0;
            end
        end
        mT = mT.*(~~L); % Eliminates small clusters
        if isempty(find(mT))
            fprintf('No sufficiently sized clusters within ROI!\n');
            stop = 1;
        end
    end
    if ~stop
    % find the max
    mT(mT==0)=NaN; %If the max threshold in the ROI is negative, we don't want to find 0 voxels.
    max_ind = find(mT==max(max(max(mT))));
    mT(~isfinite(mT))=0;
    [mx my mz] = ind2sub(t.dim,max_ind);
    max_xyz = [mx my mz];
    if config.top_voxels
        if config.top_voxels>length(find(mT)) % search space has less than the number of requested voxels
            fprintf('warning: there are less than %i voxels with T value over 0 in your search space for %s \n',config.top_voxels,names{r});
        else            
%            null = sort(unique(mT),'descend');
%Change the 0-voxels (almost entirely outside the brain, etc) into -999
%voxels so they won't get chosen.
            mT(mT==0) = -999;
            [vals,ind]=sort(mT(:),'descend');
            
%            null = sort(mT(:),'descend');
%            mT(mT<null(config.top_voxels)) = 0;
            vals(config.top_voxels+1:end) = 0; %Zero out all but the top n voxels
            mT(ind) = vals;
        end
    else
        if config.old_style
            % then you need to jump to the nearest peak
            % compute the distance to the true peaks
            peak_dist = (sum(bsxfun(@minus,peaks_xyz,max_xyz).^2,2).^(1/2));
            max_ind = peaks_ind(find(peak_dist == min(peak_dist)));
            [mx my mz] = ind2sub(t.dim,max_ind);
            max_xyz = [mx my mz];
        end
        clust_id = L(max_ind);
        mT = mT.*(L==clust_id);  %Eliminates all clusters except for the one with max peak
    end
    %limit the cluster to a sphere around the peak voxel with the configed
    %radius (unless configured not too)
    stv_inds = find(mT);
    [sx sy sz] = ind2sub(t.dim,stv_inds);
    stv_xyz = [sx sy sz];
    if ~config.whole_cluster
        stv_dist = (sum(bsxfun(@minus,stv_xyz,max_xyz).^2,2).^(1/2));
        clust_inds = stv_inds(stv_dist<=(config.roi_radius/2)); % divide by two, because we're in voxel space
    else
        clust_inds = stv_inds;
    end
    %collect the picked ROI data into variables and save ROI .mat file
    peak_T = mT(max_ind);
    mean_T = mean(mT(find(mT)));
    min_T = min(mT(find(mT)));
    n_vox = size(clust_inds,1);
    peak_xyz_vox = max_xyz;
    peak_xyz_ind = max_ind;
    peak_xyz_mm = XYZ(:,max_ind)';
    [rx ry rz] = ind2sub(t.dim,clust_inds);
    roi_XYZ = [rx ry rz];
    roi_XYZmm = XYZ(:,clust_inds)';
    roi_XYZind = clust_inds;
    if config.whole_cluster
        if config.top_voxels
            xY.descrip = 'top_voxels';
            xY.str = sprintf('top%dvoxels',config.top_voxels);
        else            
            xY.descrip = 'cluster';
            xY.str = 'cluster';
        end
    else
        xY.descrip = 'sphere';
        xY.str = sprintf('%Gmm_sphere',config.roi_radius);
    end
    xY.xyz = peak_xyz_mm;
    xY.M = t.mat;
    xY.XYZmm = roi_XYZmm';
    rname = [names{r} '_' xY.str roi_suffix];
    save(fullfile(targ,[rname '.mat']),'xY','config','peak_T','mean_T','min_T','n_vox', ...
         'peak_xyz_vox','peak_xyz_ind','peak_xyz_mm','roi_XYZ','roi_XYZmm','roi_XYZind');  
    try
    	f=fopen(fullfile(job.log_location,job.log_fname{r}),'a');
        if config.top_voxels
            fprintf(f,'%s,%i,%i,%i,%G,%G,%G,%G\n',subject,peak_xyz_mm(1),peak_xyz_mm(2),peak_xyz_mm(3),n_vox,mean_T,peak_T,min_T);
        else
            fprintf(f,'%s,%i,%i,%i,%i,%G\n',subject,peak_xyz_mm(1),peak_xyz_mm(2),peak_xyz_mm(3),n_vox,peak_T);
        end
        fclose(f);
    end
    rh = t;
    rh.fname = fullfile(targ,[rname '.img']);
    R = zeros(rh.dim);
    R(roi_XYZind) = 1;
    spm_write_vol(rh,R);
    % if you're supposed to render stuff, now is the time!
    if config.render
        render_images(study,subject,model,contrast,config,names{r},rname,filename);
    end
    end
end
% convert the .ps file, if it exists!
if exist([filename '.ps'],'file')
%    ps2pdf('psfile',[filename '.ps'],'pdffile',[filename '.pdf']);
%    system(sprintf('rm -rf %s',[filename '.ps']));
end
end


%% Render image function
function render_images(study,subj,model,con,config,roi,rname,filename)
%-------------------------------------------------------------------------%
% render_images(study,subj,model,con,roi,n2r,c2r)
% study = string, study
% subj = string, subject
% model = string, model
% con = int, contrast number from which to obtain the T
% roi = string, the ROI name
% n2r = [int, int], number of images to render simultaneously
% c2r = int, current image to render (index into the subplot)
%-------------------------------------------------------------------------%
% this script renders the ROI of interest, overlayed on the participant's
% anatomical and functional activation. The rendering is done in a subplot,
% the size of which is given by the 
%subplot(n2r(1),n2r(2),c2r);


%when picking voxels threshold is set to 1 to eliminate the p_value
%constraint. Ploting this would be really messy and uninformative.
if config.p_val_thresh == 0.5         
    thresh = 0.001;
else
    thresh = config.p_val_thresh;
end

func_trans = 0.75;
ROI_trans = 1;
color_allocation = [  1 200; ...
                    201 300; ...
                    301 301;]; % this is the color range allocated to the anatomical, functional, and ROI images respectively
                
root = char(adir(fullfile(study,subj)));
try
    a = spm_vol(char(adir(fullfile(root,'3danat','ws*.img'))));
catch
    a = spm_vol(char(adir('/mindhive/saxelab/scripts/template/rT1.nii')));
end
A = spm_read_vols(a); % anatomical
try
    f = spm_vol(char(adir(fullfile(root,'results',model,sprintf('spmT_%04i.img',con)))));
catch
    f = spm_vol(char(adir(fullfile(root,sprintf('spmT_%04i.img',con)))));
end
F = spm_read_vols(f); % functional T-map 
r = char(adir(fullfile(root,'autoROI',[rname '.mat'])));
r = load(r);
R = zeros(size(F));R(r.roi_XYZind) = 1;

% threshold the functional T-map
% => aquire degrees of freedom
null = regexp(f.descrip,'[[]]','split');
dof = str2num(null{2});
% => compute the critical T-statistic
t_crit = tinv(1-thresh,dof);
% => filter & threshold T-map

%if ~config.top_voxels
    F(~isfinite(F))=0;
    F(F<t_crit)=0;
%end


% Determine the peak T-value
MAT    = f.mat;
IMAT   = inv(MAT);
xyz = r.xY.xyz';
%xyz(4) = 1;
xyz(4,1:end) = 1;
xyz = IMAT*xyz;xyz=xyz(1:3); % xyz is the center of the roi

% produce appropriate image slices
Ax=squeeze(A(xyz(1),:,:));Fx=squeeze(F(xyz(1),:,:));Rx=squeeze(R(xyz(1),:,:));
Ay=squeeze(A(:,xyz(2),:));Fy=squeeze(F(:,xyz(2),:));Ry=squeeze(R(:,xyz(2),:));
Az=squeeze(A(:,:,xyz(3)));Fz=squeeze(F(:,:,xyz(3)));Rz=squeeze(R(:,:,xyz(3)));
Ax=rot90(Ax);Fx=rot90(Fx);Rx=rot90(Rx);
Ay=rot90(Ay);Fy=rot90(Fy);Ry=rot90(Ry);

Ay=flipdim(Ay,2);Fy=flipdim(Fy,2);Ry=flipdim(Ry,2);
Az=flipdim(Az,1);Fz=flipdim(Fz,1);Rz=flipdim(Rz,1);
% compute alpha mapping
Axa = ones(size(Ax));Aya = ones(size(Ay));Aza = ones(size(Az));
Fxa = zeros(size(Fx));Fxa(find(Fx)) = func_trans;
Fya = zeros(size(Fy));Fya(find(Fy)) = func_trans;
Fza = zeros(size(Fz));Fza(find(Fz)) = func_trans;
Rxa = zeros(size(Rx));Rxa(find(Rx)) = ROI_trans;
Rya = zeros(size(Ry));Rya(find(Ry)) = ROI_trans;
Rza = zeros(size(Rz));Rza(find(Rz)) = ROI_trans;

% now, store the images in a cell array to make indexing easier
% imgs{a,b,c} >
%       a, image (anat, func, roi)
%       b, axis (x, y, z)
%       c, type (image, alpha map)

% oh right, and the y-images have to be flipped Left to Right while the z
% images have to be flipped up > down

imgs{1,1,1} = Ax;imgs{1,1,2} = Axa;
imgs{1,2,1} = Ay;imgs{1,2,2} = Aya;
imgs{1,3,1} = Az;imgs{1,3,2} = Aza;
imgs{2,1,1} = Fx;imgs{2,1,2} = Fxa;
imgs{2,2,1} = Fy;imgs{2,2,2} = Fya;
imgs{2,3,1} = Fz;imgs{2,3,2} = Fza;
imgs{3,1,1} = Rx;imgs{3,1,2} = Rxa;
imgs{3,2,1} = Ry;imgs{3,2,2} = Rya;
imgs{3,3,1} = Rz;imgs{3,3,2} = Rza;

% oh right, and the y-images have to be flipped Left to Right while the z
% images have to be flipped up > down


% now we need to generate a colormap; to force it to use multiple
% colormaps, we have to catenate a bunch of colormaps together
% the anatomical will be displayed in GREY
% ... guh...
% actually, this isn't that hard. Let's allocate the range 0-200 to the
% anatomical, 201-300 to the functional, and 301 to the ROI!
% so let's set all the images to their appropriate ranges
ca=color_allocation; %to make it easier to type
for axis = 1:3
    for img = 1:3
        n=reshape(imgs{img,axis,1},[],1);n=n+(-1*min(n));n=n/max(n);
        n=n*(ca(img,2)-ca(img,1));
        n=n+ca(img,1);
        imgs{img,axis,1}=round(reshape(n,size(imgs{img,axis,1})));
    end
end
CMap = [gray(length(ca(1,1):ca(1,2))); hot(length(ca(2,1):ca(2,2))); 0 0 1];
colormap(CMap);
% let's not do the subplot thing; it's useless. why not catenate the images
% together instead
%dims = a.dim;
imsz = [size(Ax,1)+size(Az,1) size(Ax,2)+size(Ay,2)];
im{1,1} = zeros(imsz);im{1,2} = zeros(imsz);
im{2,1} = zeros(imsz);im{2,2} = zeros(imsz);
im{3,1} = zeros(imsz);im{3,2} = zeros(imsz);
col_inds = [1 size(imgs{1,1,1},2); size(imgs{1,1,1},2)+1 size(imgs{1,1,1},2)+size(imgs{1,2,1},2); 1 size(imgs{1,3,1},2)];
row_inds = [1 size(imgs{1,1,1},1); 1 size(imgs{1,2,1},1); size(imgs{1,1,1},1)+1 size(imgs{1,1,1},1)+size(imgs{1,3,1},1)];
for imz = 1:size(imgs,1)
    for ax = 1:size(imgs,2)
        im{imz,1}(row_inds(ax,1):row_inds(ax,2),col_inds(ax,1):col_inds(ax,2)) = imgs{imz,ax,1};
        im{imz,2}(row_inds(ax,1):row_inds(ax,2),col_inds(ax,1):col_inds(ax,2)) = imgs{imz,ax,2};
    end
end
h(1) = image(im{1,1});
hold on
h(2) = image(im{2,1},'AlphaData',im{2,2});
h(3) = image(im{3,1},'AlphaData',im{3,2});
set(h(1),'CDataMapping','direct');
set(h(2),'CDataMapping','direct');
set(h(3),'CDataMapping','direct');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

%-------------------------------------------------------------------------%
% LABELING
% let's draw lines? 
% it seems that the horizontal line is given by the z-coordinate, or rather
% its compliment
line([1 size(im{1,1},2)],[a.dim(3)-xyz(3) a.dim(3)-xyz(3)],'Color','r');
% the veritcal line is merely the y-coordinate
line([xyz(2) xyz(2)],[1 size(im{1,1},1)],'Color','r');
hold off
s2p = sprintf('%s\n%s\n%s\nSize: %i\nPeak T: %.3f\nCoords: %i %i %i\n',subj,model,roi,size(r.xY.XYZmm,2),F(xyz(1),xyz(2),xyz(3)),r.xY.xyz);
s2p=strrep(s2p,'_','\_');
text(size(Ax,2)+10,round((size(im{1,1},1))*2/3),s2p,'Color','w');
drawnow;
print(gcf,filename,'-r150','-append','-dpsc2')
end
