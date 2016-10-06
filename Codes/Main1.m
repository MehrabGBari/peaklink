clc
clear all

config

ExperimentID = 'exp01';
mzXMLfileName1 = '20110922_EXQ4_NaNa_SA_YeastEasy_Labelfree_01.mzXML';
mzXMLfileName2 = '022008_F4_DTree_1.mzXML';
mzXMLfileName3 = '022008_F5_DTree_1.mzXML';

maxquantFolder1 = '20110922_EXQ4_NaNa_SA_YeastEasy_Labelfree_01';
maxquantFolder2 = '022008_F4_DTree_1';
maxquantFolder3 = '022008_F5_DTree_1';

getPath
tic
preProcessing
toc
runNotification('Processing Done!');

tic
[Parameters, PP_Parameter, TrainingId_AB, TrainingId_AC, ~] = generateParameters(INTERVAL_Matrix, Information_Matrix_Total, retentiont_A, retentiont_B, retentiont_C, retentiontl1_A, retentiontl1_B, retentiontl1_C, XICs_A, XICs_B, XICs_C, IntervalListBig);
TestModel = buildATARAKLModel(Information_Matrix_Total, INTERVAL_Matrix, Parameters, IntervalListBig, retentiontl1_A, retentiontl1_B, retentiontl1_C, retentiont_A, retentiont_B, retentiont_C, XICs_A, XICs_B, XICs_C, PP_Parameter);
runNotification('Model Building Done!');

 [scoreErrRate, nIntersectErrRate, nTestingErrRate, pos01, posErr] = alignResEvaErrRate(TrainingId_AB, TrainingId_AC, TestModel, Information_Matrix_Total, INTERVAL_Matrix);
 [scoreDetRate, nIntersectDetRate, nTestingDetRate, pos02] = alignResEvaDetRate(TrainingId_AB, TrainingId_AC, TestModel, Information_Matrix_Total, INTERVAL_Matrix);
 

 save([save_path,ExperimentID, '_',  'ErroRate'], 'scoreErrRate', 'nIntersectErrRate', 'nTestingErrRate', 'pos01', 'posErr')
 save([save_path,ExperimentID, '_', 'DetReat'], 'scoreDetRate', 'nIntersectDetRate', 'nTestingDetRate', 'pos02')

timeConsumed = toc
logging
runNotification('Done!');
toc
