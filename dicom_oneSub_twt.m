function study = dicom_oneSub_twt(study)

tic
% establish existing variables
studyDir    = '/Users/toddt/Dropbox/SaxelabProcessing/fMRI_DataAnalysisTraining';
subjects    = {'SAX_DAT_01'};

% when doing this for all of the subjects, the for loop will need to start
% here

% create directories to store the converted data and corresponding variables 
subj        = subjects{1};
subjDir     = makeSubjDirs(studyDir,subj);
dicomDir    = subjDir{1};
resultsDir  = subjDir{2};
roiDir      = subjDir{3};
anatDir     = subjDir{4};
origDir     = subjDir{5};
scoutDir    = subjDir{6};
boldDir     = subjDir{7};
repDir      = subjDir{8};

% find all of the dicoms and create a run information report
[serDes, seqName, run] = findDicoms12(dicomDir, repDir);

% determine which dicoms are which runs (e.g. functional, anatomical, etc.)
cd(dicomDir) % make sure working from the dicom directory

% find the functional runs
mocoRuns = run(strcmp('MoCoSeries',serDes));
funcRuns = [];
if ~isempty(mocoRuns)
    funcRuns = mocoRuns - 1; % the functional run is always the run before the motion corrected run
else
    funcRuns = run(strcmp('epfid2d1_64',deblank(seqName)));
    if isempty(funcRuns)
        funcRuns = run(strcmp('epfid2d1',deblank(seqName)));
    end
end

% find the anatomical runs 
anatRuns = run(strncmp('T1_MPRAGE',serDes,9));
if isempty(anatRuns)
    anatRuns = run(strncmp('T1_MPRAG',serDes,8));
    if isempty(anatRuns)
        anatRuns = run(strcmp('MPRAG',serDes,5));
    else
        anatRuns = []; % if there are no anatomical runs, store an empty vector in the variable 
    end
end

% find everything else ('scoutRun' is anything that isn't moco, func, or anat)
scoutRuns = setdiff(run, [mocoRuns funcRuns anatRuns]);

% convert all of the runs and move them to the correct directories

% anatomicals first - It is possible there will not be anatomical
% runs. This script assumes that the user will know when there are and are
% not anatomicals in the dicoms and does not generate an error message.
if ~isempty(anatRuns)
    convertDicom(anatRuns(end));
    system(sprintf('cp %s %s',[dicomDir sprintf('/*-%04i-*.[^d]*',anatRuns(end))],origDir));
    system(sprintf('mv %s %s',[dicomDir sprintf('/*-%04i-*.[^d]*',anatRuns(end))],anatDir));
end

% functionals
convertDicom(funcRuns);
for i = funcRuns
    cbd = fullfile(boldDir, sprintf('%03i',i)); mkdir(cbd); % current bold directory
    system(sprintf('mv %s %s',[dicomDir sprintf('/*-%04i-*.[^d]*',i)],cbd));
end

% all other runs
convertDicom(scoutRuns);
for j = scoutRuns
    csd = fullfile(scoutDir,sprintf('%03i',j)); mkdir(csd); % current scout directory
    try
        system(sprintf('mv %s %s', [dicomDir sprintf('/*-%04i-*.[^d]*',j)],csd));
    end
end

cd(studyDir)

% now the thing I would like to add to this script is to re-create the
% behavioural file somehow. Count the TRs and see if we can match up the
% bold directories to the correct behavioural file and save that
% information in the structure. 
toc

end
