for i = 8:15

study.name='fMRI_DataAnalysisTraining';
study.path='/Users/toddt/fMRI_DataAnalysisTraining';
study.subjects={['SAX_DAT_' num2str(i,'%02d')]};

%study = dicom_oneSub_twt(study)
%study = prep_oneSub_2017(study)
%saxelab_generate_art_2017(study.path,study.subjects{1},[],{'overwrite',1})
%artreport = saxelab_art_mot_report_2017(study.path,study.subjects{1})
saxelab_compcor_2017(study.path,study.subjects{1},[],5,1,'')
%saxelab_model_2017(study.path,study.subjects{1},'tom',[09 10 11 12], {'add_art',1, 'add_cc', 1})

switch i
case 6
bolds = [11 12 13 14]
case {9,11,15}
bolds = [10 11 12 13]
case 8
bolds = [9 11 12 13]
otherwise
bolds = [9 10 11 12]
end

saxelab_model_2017(study.path,study.subjects{1},'tom',bolds, {'add_art',1, 'add_cc',1})
end
