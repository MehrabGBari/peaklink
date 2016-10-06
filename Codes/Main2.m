close all;clear all;clc

%[Parameters, PP_Parameter, TrainingId_AB, TrainingId_AC, ~] = generateParameters(INTERVAL_Matrix, Information_Matrix_Total, retentiont_A, retentiont_B, retentiont_C, retentiontl1_A, retentiontl1_B, retentiontl1_C, XICs_A, XICs_B, XICs_C, IntervalListBig);

%% AB

% Train Model sellection
 
TrainingId_ABTmp = genTrainingId(INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,5)/2);
k=1;
tic

m=floor(0.435*length(TrainingId_ABTmp));
TrainingId_AB = TrainingId_ABTmp(1:m);
TrainingInformation_AB = trainModel(retentiontl1_A, retentiont_A, XICs_A, retentiontl1_B, retentiont_B, XICs_B, Information_Matrix_Total, IntervalListBig, INTERVAL_Matrix, TrainingId_AB, 1);

%AT analysis
PP_Parameter.PP_AB = genPP_parameter(retentiont_A, retentiont_B, Information_Matrix_Total, TrainingId_AB, 1);
AT_Model_Training_AB = Generate_Training_ATScores(TrainingInformation_AB, PP_Parameter, 1);
[Parameters,Threshold] = BuildATModelParameters(TrainingInformation_AB,AT_Model_Training_AB);
%AW Analysis

[AW_Model_Training_AB,Trulable]=Generate_Training_AWScores(TrainingInformation_AB,PP_Parameter, 1,Parameters,Threshold);%taks some times

[Parameters,PAWScores,SVMStruct] = BuildAWModelParameters_1(TrainingInformation_AB,AW_Model_Training_AB,Parameters);


%Calculating AT score
TestingId_AB = TrainingId_ABTmp(m+1:end);
TestingInformation_AB = TestModel(retentiontl1_A, retentiont_A, XICs_A, retentiontl1_B, retentiont_B, XICs_B, Information_Matrix_Total, IntervalListBig, INTERVAL_Matrix, TestingId_AB, 1);
Accuracy_withAT = Calculate_Testing_ATScores(TestingInformation_AB, PP_Parameter, 1);


[AW_Model_Testing_AB,Ist_step]=Generate_Training_AWScores_test(TestingInformation_AB,PP_Parameter, 1,Parameters,Threshold);%taks some times

[TestAccuracy_SVM] = BuildAWModelParameters_1_test(TestingInformation_AB,AW_Model_Testing_AB,Parameters,SVMStruct);

toc

%AW Analysis
tic
[AW_Model_Training_AB,Trulable]=Generate_Training_AWScores(TrainingInformation_AB,PP_Parameter, 1,Parameters,Threshold);%taks some times
toc
[Parameters,PAWScores,SVMStruct] = BuildAWModelParameters_1(TrainingInformation_AB,AW_Model_Training_AB,Parameters);

%%
%%Testing
TestingId_AB = TrainingId_ABTmp(m+1:end);
TestingInformation_AB = TestModel(retentiontl1_A, retentiont_A, XICs_A, retentiontl1_B, retentiont_B, XICs_B, Information_Matrix_Total, IntervalListBig, INTERVAL_Matrix, TestingId_AB, 1);

tic
[AW_Model_Testing_AB,Trulable]=Generate_Training_AWScores_test(TestingInformation_AB,PP_Parameter, 1,Parameters,Threshold);%taks some times
toc
[PAWScores] = BuildAWModelParameters_1_test(TestingInformation_AB,AW_Model_Testing_AB,Parameters,SVMStruct);

