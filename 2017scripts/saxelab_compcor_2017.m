function saxelab_compcor_2017(study_dir,subjects,bolds,pcnum,use_art,art_prefix)
%This script runs PCA on signal extracted from white matter/csf masks, and
%saves the top [pcnum] components as nuisance regressors, to be included in
%subsequent modeling.

%%%INPUTS
% study: A string w/ the study name (e.g. 'KMVPA')
% subjects: a single string/cell array w/ subject strings (e.g. 'SAX_KMVPA_02', or {'SAX_KMVPA_02','SAX_KMVPA_03'})
% bolds: specify which bold dirs to include (e.g. [11 13 15 17], or {[11 13 15 17],[13 15 19 21]}
% pcnum: number of principle components to save as nuisance regressors (5 is common).
% art_tps: 0 if you don't want to acct for artifact tps (IDed w/ art), 1 if
% you do
% art_prefix: art_prefix for art_reg files (string)- note, you need this
% regardless of whether you're accounting for arts in compcor
% calculations, bc it will create a file w/ your art + cc regressors. 
% Example command (1 subject): saxelab_compcor('KMVPA','SAX_KMVPA_02',[11 13 15 17 19],5,1,'2mm3std')
% Example command (multiple subjects): saxelab_compcor('KMVPA',{'SAX_KMVPA_02','SAX_KMVPA_03'},{[11 13 15 17 19],[13 17 19 21 23]},5,0,'')

%%%OUTPUTS
%In {subject}/compcor:
% ccprep.mat file: w/ art regressor info and numRuns
% ccresults.mat file: w/ mean vol, noisepool, ccregs, and pcnum variables
% figuredir w/ picture of mean vol and noise pool
% in each bolddir: cc regressors art file + .txt file


%% Step 1: Setup
addpath('/mindhive/saxelab/scripts/')
addpath('/mindhive/saxelab/scripts/GLMdenoise-1.4/utilities/')

%Study directory
%study_dir = adir(['/mindhive/saxelab*/' study]);
%if ~iscell(study_dir)
if ~isdir(study_dir)
    error('Could not locate study. Please provide full path');
end
%study_dir = study_dir{1}; %have to do this in order to cd
cd(study_dir);

% Check that the subject <subjects> is type CHAR or CELL
%if ~iscell(subjects)
    if ischar(subjects)
        subjects = {subjects};
    else
%        error('subjects must be a STRING or CELL array of STRINGs');
        error('Subject must be a STRING');
    end
%end

% Check that the bold directories <bolds> are type DOUBLE or CELL
% if ~isnumeric(bolds)&&~iscell(bolds)
%     error('BOLDS must be an array of DOUBLES or CELL array of arrays of DOUBLEs');
%end

if ~isnumeric(bolds)
    error('BOLDS must be an array of DOUBLES');
end

if isnumeric(bolds)
    bolds = {bolds};
end

% if iscell(bolds)
%     if length(bolds)==1
%         bolds = repmat(bolds,1,length(subjects));
%     elseif length(bolds)~=length(subjects)
%         error('Number of sets of bold directories specified does not match the number of subjects');
%     end
% end

%% Serialize the processing of subjects
% if length(subjects)>1
%     for s = 1:length(subjects)
%         saxelab_compcor(study,subjects{s},bolds{s},pcnum,use_art,art_prefix)
%     end
%     return
% end

%% Iterate through subjects, pulling out both design matrix info and fMRI data
subject  = subjects{1};
bolds = bolds{1};

fprintf(['Processing data for ' subject '\n'])

%subject_dir = adir(['/mindhive/saxelab*/' study '/' subject]);
%subject_dir = subject_dir{1}; %have to do this in order to cd
subject_dir = fullfile(study_dir, subject);
cd(subject_dir);

if ~exist('compcor','dir') == 1
    mkdir('compcor');
end

% make sure bolds are in the correct orientation, just in case
bolds = reshape(bolds,1,[]);
% check that you can find the bold directories
bdirs = {};

if isempty(bolds)
    %        bdirs = adir(fullfile(subjects{s},'bold','0*'));
    bdirs = cellstr(spm_select('FPList',fullfile(study_dir,subjects,'bold'),'dir'));
    %        if ~iscell(bdirs)
    if isempty(bdirs{1})
        %             fprintf('could not locate any bold directories for %s, skipping subject.\n',subjects{s});
        %             continue
        error('Could not locate any bold directories for subject.');
    end
else
    
    for i = bolds
        %    if ~iscell(adir(fullfile(subject_dir,'bold',sprintf('%03d',i))))
        if ~isdir(fullfile(subject_dir,'bold',sprintf('%03d',i)))
            error(sprintf('Could not locate bold directory %03d',i));
        end
        %    bdirs(end+1) = adir(fullfile(subject_dir,'bold',sprintf('%03d',i)));
        bdirs{end+1} = fullfile(subject_dir,'bold',sprintf('%03d',i));
    end
end


numRuns = length(bdirs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 1: Create data variable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Creating data variable for ' subject '\n'])
    cd(study_dir);
    cd(subject);

    data = cell(1,numRuns);
    use_tps = cell(1,numRuns);
    
    if isempty(art_prefix)
        art_prefix = '';
    else
        art_prefix = [art_prefix '_'];
    end
    
    if use_art == 1
        art_reg = cell(1,numRuns);
        art_tps = cell(1,numRuns);
    end

    for bdir = 1:length(bdirs)
        cd(bdirs{bdir});

%        vols = adir('swrf*.img');
        vols = cellstr(spm_select('FPList',bdirs{bdir},'^swrf.*.img'));
        all_tps = (1:length(vols))';
        
        if use_art == 1
%            art_reg_info = adir([art_prefix 'art_regression_outliers_swrf*.mat']); 
            art_reg_info = cellstr(spm_select('FPList',bdirs{bdir}, ['^' art_prefix 'art_regression_outliers_swrf.*.mat'])); 
            load(art_reg_info{1});
            art_reg{bdir} = R;
        
            use_vols = length(vols)-size(R,2); %how many vols
            use_tps{bdir} = zeros(use_vols,1); %which tps to use
            art_tps{bdir} = zeros(size(art_reg{bdir},2),1);
              
            if size(art_reg{bdir},2) > 0

                for h = 1:size(art_reg{bdir},2)
                    art_tps{bdir}(h,1) = find(art_reg{bdir}(:,h));
                end

                use_tps{bdir} = setdiff(all_tps,art_tps{bdir});

                if size(use_tps{bdir}) ~= use_vols
                    error('number of TPs must be equal to number of images minus number of arts');
                end
            else
                use_tps{bdir} = all_tps;
            end
        
        else
            use_vols = length(vols);
            use_tps{bdir} = all_tps;
        end
        
        data{bdir} = zeros(79,95,79,use_vols); %should match ips, but if scanner stopped short.. minus num artTPs
        
        for vol_num = 1:size(use_tps{bdir},1) %for every image, read in brain data
            vol = use_tps{bdir}(vol_num);
            v=spm_vol(vols{vol}); %read in one image
            d=spm_read_vols(v);  %read in one image
            size(d); %should be 79x95x79 (size of one image)

            %save d as xyz for vol; vol is 4th (time) dimension
            data{bdir}(:,:,:,vol_num) = d;
            data{bdir} = single(data{bdir}); %example dataset is single so we'll do that too.
        end
            
    end

    cd([subject_dir '/compcor']);
    
    if use_art == 1
        save([subject '.ccprep.mat'],'data','art_reg','use_tps','art_tps','numRuns');
    else
        save([subject '.ccprep.mat'],'data','use_tps','numRuns');
    end
    clear data;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 2: Create noise pool (eroded white matter mask + csf)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Creating noisepool for ' subject '\n'])
    cd([subject_dir '/mask']);

    %load WM mask
    %WM = spm_vol('white_matter_mask.img');
    WM = spm_vol('white_matter_mask.nii');
    wmmask_full = spm_read_vols(WM);
    wmmask_full = logical(wmmask_full); %3D binary image matrix
    
    %initialize header info for noisepool image
%    np = adir('ws0-*.img');
    np=cellstr(spm_select('FPList',fullfile(subject_dir, '3danat'), '^ws0-.*.img'));
    np = spm_vol(np{1});
    np.fname = 'noisepool.img';

    %Erode WM map by two voxels in each direction, to minimize partial
    %voluming w/ grey matter
    if ~(isnumeric(wmmask_full) || islogical(wmmask_full)) || numel(size(wmmask_full))~=3
        disp('Input is not a 3D image matrix.');
        return;
    end;

    if numel(wmmask_full)~=numel(wmmask_full(ismember(wmmask_full,[0 1])))
        disp('Input image must be a binary mask.');
        return;
    end;

    dims = size(wmmask_full);
    wmmask = wmmask_full;

    for i = 2:(dims(1)-1)
        for j = 2:(dims(2)-1)
            for k = 2:(dims(3)-1)
                if sum(sum(sum(wmmask_full((i-1):(i+1),(j-1):(j+1),(k-1):(k+1)))))~=27
                    wmmask(i,j,k)=0;
                end;
            end;
        end;
    end;

    %if it exists, load csf mask
%    if exist('csf_mask.img','file') > 0
%        csf = spm_vol('csf_mask.img');
    if exist('csf_mask.nii','file') > 0
        csf = spm_vol('csf_mask.nii');
        csfmask_full = spm_read_vols(csf);
        csfmask_full = logical(csfmask_full); %3D binary image matrix
        
        dims2 = size(csfmask_full);
        csfmask = csfmask_full;

        for l = 2:(dims2(1)-1)
            for m = 2:(dims2(2)-1)
                for n = 2:(dims2(3)-1)
                    if sum(sum(sum(csfmask_full((l-1):(l+1),(m-1):(m+1),(n-1):(n+1)))))~=27
                        csfmask(l,m,n)=0;
                    end;
                end;
            end;
        end;

        %take conjunction of eroded wm mask and csf mask as noise pool
        noisepool = wmmask | csfmask ;

        cd([subject_dir '/compcor']);
        save('wmmask','wmmask');
        save('csfmask','csfmask');
        save('noisepool','noisepool');
        nopo = spm_write_vol(np,noisepool);
        clear wmmask csfmask noisepool;
    else
        noisepool = wmmask; %otherwise, just use the wmmask

        cd([subject_dir '/compcor']);
        save('wmmask','wmmask');
        save('noisepool','noisepool');
        nopo = spm_write_vol(np,noisepool);
        clear wmmask csfmask noisepool;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Step 3: Run CompCor!
    %%%%%%%%%%%%%%%%%%%%%%%%
    %NOTE: subscript.m, projectionmatrix.m, unitlengthfast.m, and svds.m
    % live in saxelab/scripts/GLMdenoise-1.4/utilities)
    fprintf(['Running compcor on subject ' subject '\n'])
    cd([subject_dir '/compcor']);

    %Load in the reformatted data
    load([subject '.ccprep.mat']) %this will load numRuns variable + art_reg variable
    load('noisepool.mat');

    % calc
    tr = 2;
    numruns = length(bdirs);
    dataclass = class(data{1});  % will always be 'single'
    is3d = size(data{1},4) > 1;
    if is3d
      dimdata = 3;
      dimtime = 4;
      xyzsize = sizefull(data{1},3);
    else
      dimdata = 1;
      dimtime = 2;
      xyzsize = size(data{1},1);
    end
    numvoxels = prod(xyzsize);

    maxpolydeg = zeros(1,numruns);
    for p=1:numruns
      maxpolydeg(p) = round(((size(data{p},dimtime)*tr)/60)/2);
    end

    volcnt = cellfun(@(x) size(x,dimtime),data);
    meanvol = reshape(catcell(2,cellfun(@(x) squish(mean(x,dimtime),dimdata),data,'UniformOutput',0)) ...
        * (volcnt' / sum(volcnt)),[xyzsize 1]);

    % determine noise regressors
    fprintf('*** glmdenoise: calculating noise regressors. ***\n');
    ccregressors = {};
    for p=1:length(data)

        % extract the time-series data for the noise pool 
        temp = subscript(squish(data{p},dimdata),{noisepool ':'})';  % time x voxels

        % project out polynomials from the data
        temp = projectionmatrix(constructpolynomialmatrix(size(temp,1),0:maxpolydeg(p))) * temp;

        % unit-length normalize each time-series (ignoring any time-series that are all 0)
        [temp,len] = unitlengthfast(temp,1);
        temp = temp(:,len~=0);

        % perform SVD and select the top PCs
        [u,s,v] = svds(double(temp*temp'),pcnum);
        u = bsxfun(@rdivide,u,std(u,[],1));  % scale so that std is 1
        ccregressors{p} = cast(u,dataclass);

    end
    clear temp len u s v;

    % the number of PCs is specified by the user
    fprintf('*** user-specified number of PCs is %d. ***\n',pcnum);

    ccresults.meanvol = meanvol;
    ccresults.noisepool = noisepool;
    ccresults.ccregressors = ccregressors;
    ccresults.pcnum = pcnum;

    cd([subject_dir '/compcor']);
    save([subject '.ccresults.mat'], 'ccresults');

    %MAKE FIGS
%     ccresults.figuredir = 'compcor_figures';
%     mkdir(ccresults.figuredir);
%     % write out image showing mean volume (of first run)
%     imwrite(uint8(255*makeimagestack(ccresults.meanvol,1)),gray(256),fullfile(ccresults.figuredir,'MeanVolume.png'));
% 
%     % write out image showing noise pool
%     imwrite(uint8(255*makeimagestack(ccresults.noisepool,[0 1])),gray(256),fullfile(ccresults.figuredir,'NoisePool.png'));

    clear ccresults meanvol noisepool ccregressors pcnum data
    
    %%%%%%%%%%%%%%%%%%%%%
    %Step 4: Save CC regs
    %%%%%%%%%%%%%%%%%%%%%
    fprintf(['Saving CC regs as regressors for subject ' subject '\n'])
    cd([subject_dir '/compcor/']);

    %delete data variable, to save space
    load([subject '.ccprep.mat']);
    if exist('data','var')
        if use_art == 1
            save([subject '.ccprep.mat'],'art_reg','use_tps','art_tps','numRuns');
            clear art_reg
        else
            save([subject '.ccprep.mat'],'use_tps','numRuns');
        end
    end

    load([subject '.ccresults.mat']);
    
    %%linearly interpolate surrounding timepoints to fill in missing TPs.
    if use_art == 1
        for b = 1:size(ccresults.ccregressors,2) %for each bolddir
            p = ccresults.ccregressors{b};

            a = art_tps{b};

            if size(a,1) > 0
                for i=1:length(a) %for each art TP
                    p = [p(1:a(i)-1,:); zeros(1,ccresults.pcnum); p(a(i):end,:)]; %find artTP row, fill with zero (temporary)
                end
                ccresults.ccregressors{b} = p;

                leftover = find(all(p==0,2)); %all chunks. but chunks aren't all next to eachother.
                D = ~~([1 diff(leftover' - (1:length(leftover')))]); %find where chunks start
                starts = find(D); %index starts
                spots = zeros(1,length(starts));

                i = 1; 
                while i<length(starts)
                    spots(i) = 1;
                    if starts(i) < length(D)-1
                        counter = 1;
                        while D(starts(i)+counter) == 0 
                            spots(i) = spots(i) + 1;
                            counter = counter + 1;
                        end
                    elseif starts(i) == length(D)-1;
                        if D(starts(i) + 1) == 0 
                            spots(i) = spots(i) + 1;
                        end
                    end
                    i=i+1;
                end
                    
                if i==length(starts) %if last value in starts
                    spots(i) = length(D) - starts(i) + 1;
                end

                for i=1:length(starts)
                    for reg = 1:ccresults.pcnum
                        %set beg tp (tp before art, or 0)
                        if a(starts(i)) == 1 %if first tp is artreg
                            beg = p(a(starts(i)),reg); %set to 0
                        else
                            beg = p(a(starts(i))-1,reg); %set to prev tp
                        end
                        
                        %set last tp (tp after art, or 0)
                        if a(starts(i))+spots(i) == all_tps(end)+1
                            last = p(a(starts(i))+spots(i)-1,reg); %set to 0
                        else
                            last = p(a(starts(i))+spots(i),reg); %set to next tp
                        end
                        
                        %interpolate b/w beg and last tps
                        if a(starts(i)) == 1
                            p(starts(i):spots(i)+1,reg) = linspace(beg,last,spots(i)+1)';
                        elseif a(starts(i))+spots(i) == all_tps(end)+1
                            p(a(starts(i))-1:a(starts(i))+spots(i)-1,reg) = linspace(beg,last,spots(i)+1)';
                        else
                            p(a(starts(i))-1:a(starts(i))+spots(i),reg) = linspace(beg,last,spots(i)+2)';
                        end
                    end                 
                end

                ccresults.ccregressors{b} = p;
            end

        end %cycle through bdirs
        save([subject '.ccresults.mat'], 'ccresults');
    end





%     % make sure bolds are in the correct orientation, just in case
%     bolds = reshape(bolds,1,[]);
%     % check that you can find the bold directories
%     bdirs = {};
%     for i = bolds
% %        if ~iscell(adir(fullfile(subject_dir,'bold',sprintf('%03d',i))))
%         if ~isdir(fullfile(subject_dir,'bold',sprintf('%03d',i)))
%             error(sprintf('Could not locate bold directory %03d',i));
%         end
% %        bdirs(end+1) = adir(fullfile(subject_dir,'bold',sprintf('%03d',i)));
%         bdirs{end+1} = fullfile(subject_dir,'bold',sprintf('%03d',i));
%     end

    cc_reg = cell(1,length(bdirs));

    for bdir = 1:length(bdirs)
        cd(bdirs{bdir});

        if use_art == 1
            %first columns of cc_regs file will be artifact regs (from gen_art)
%            art_reg_info = adir([art_prefix 'art_regression_outliers_swrf*.mat']);
            art_reg_info = cellstr(spm_select('FPList',bdirs{bdir}, ['^' art_prefix 'art_regression_outliers_swrf.*.mat']));
            load(art_reg_info{1});
            cc_reg{bdir} = R;
        end

        %last 5 columns of cc_regs file will be CC regs (from GLMDenoise)
        pc_num = ccresults.pcnum;
        cc_regressors = ccresults.ccregressors{bdir}(:,1:pc_num);

        clear 'R';
        R = [cc_reg{bdir} cc_regressors]; %art outliers, 5 pcs

        bolddir = regexp(pwd,'/','split');
        bolddir = bolddir{7};

        %art and cc regs file
        save([art_prefix 'art_regression_outliers_and_' num2str(pc_num) 'ccregs_swrf_0' bolddir '-00001-000001-01.mat'],'R')

        %cc regs file
        dlmwrite([art_prefix num2str(pc_num) 'cc_f0-0' bolddir '-00001-000001-01.txt'],cc_regressors,'precision',5,'delimiter','\t')

    end %loop through bdirs

end