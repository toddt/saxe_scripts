study.name='fMRI_DataAnalysisTraining';
study.path='/Users/toddt/fMRI_DataAnalysisTraining';
study.subjects={'SAX_DAT_01'};

%study = dicom_oneSub_twt(study)
study = prep_oneSub_twt(study)