function rename_arts_2017(study,subjects,bolds,custom_prefix,overwrite)


if nargin<1
    help rename_arts_2017;
    return;
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
    error('Could not find study. Please use full path.');
end


if nargin==1||isempty(subjects)
    % find all subjects which match.
    %search_string = fullfile(study,'SAX*');
    %    subjects = adir(search_string);
    subjects = cellstr(spm_select('FPList',study,'dir','^SAX.*'));
%    if ~iscell(subjects)
    if isempty(subjects{1})
        error('Could not locate any subjects.');
    end
    fprintf('Found %i subjects to analyze.\n',length(subjects));
% elseif ischar(subjects) && any(strfind(subjects,'*'))
%     % a filter was given, find all subjects which match.
%     search_string = fullfile(study,subjects);
%     subjects = adir(search_string);
%     if ~iscell(subjects)
%         error('Could not locate any subjects.');
%     end
%     fprintf('Found %i subjects to analyze.\n',length(subjects));
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

% Check that the bold directories if defined <bolds> are type DOUBLE or CELL
if ~exist('bolds','var')
    bolds = [];
end
if ~isempty(bolds)
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
end

for s = 1:length(subjects)
    % look for bold folders
    if isempty(bolds)
%         bdirs = adir(fullfile(subjects{s},'bold','0*'));
%         if ~iscell(bdirs)
%             fprintf('could not locate any bold directories for %s, skipping subject.\n',subjects{s});
%             continue
%         end
        bdirs = cellstr(spm_select('FPList',fullfile(subjects{s},'bold'), 'dir', '^0.*'));
        if isempty(bdirs{1})
            error(['Could not locate any bold directories for subject: ' subjects{s}]);
        end
        
    else
        bdirs = {};
        for i = bolds{s}
%             if iscell(adir(fullfile(subjects{s},'bold',sprintf('%03d',i))))
%                 bdirs(end+1) = adir(fullfile(subjects{s},'bold',sprintf('%03d',i)));
%             end           
            tmp = cellstr(spm_select('FPList',fullfile(subjects{s},'bold'),'dir',sprintf('^%03d',i)));
            if isempty(tmp{1})
                error(['Could not locate specified bold directories for subject: ' subjects{s}]);
            else
                bdirs(end + 1) = tmp;
            end
         end
%         if length(bdirs) < length(bolds{s})
%             fprintf('could not locate specified bold directories for %s, skipping subject.\n',subjects{s});
%             continue
%         end
    end
    for i = 1:length(bdirs)
%         artfiles = adir(fullfile(bdirs{i},'art_regression_outliers*.mat'));
%         if iscell(artfiles)
          artfiles = cellstr(spm_select('FPList', bdirs{i}, '^art_regression_outliers.*mat'));
          if ~isempty(artfiles{1})
            newnames = cellfun(@(x) [bdirs{i} '/' custom_prefix '_' x], cellfun(@(y) y{end}, regexp(artfiles,'/','split'),'UniformOutput',0),'UniformOutput',0);
            if exist(newnames{1},'file')&~overwrite
                tmp = questdlg('Old custom name art files found! Delete?\n this will apply to all matching files','art files found','Yes','No','No');
                if ~strcmpi(tmp,'Yes');
                    return;
                else
                    overwrite = 1;
                end
            end
            for j = 1:length(artfiles)
                system(['mv -v ' artfiles{j} ' ' newnames{j}]);
            end
        end
    end
end