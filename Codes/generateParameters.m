
function [Parameters, PP_Parameter, TrainingId_ABTmp, TrainingId_ACTmp, TrainingId_BCTmp] = generateParameters(INTERVAL_Matrix, Information_Matrix_Total, retentiont_A, retentiont_B, retentiont_C, retentiontl1_A, retentiontl1_B, retentiontl1_C, XICs_A, XICs_B, XICs_C, IntervalListBig )


% Diagnostic Code Beging ----------------------------------------------------
% 
% IDMwithMS2=[INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,5)/2, INTERVAL_Matrix(:,9)/3];
% IDVwithMS2=sum(IDMwithMS2,2);
% TrainingId_train=find(IDVwithMS2==3);

% Information_Matrix_Total = Information_Matrix_Total(TrainingId_train, :);
% XICs_A = XICs_A(:, (TrainingId_train-1)*8+1:(TrainingId_train-1)*8+8);
% XICs_B = XICs_B(:, (TrainingId_train-1)*8+1:(TrainingId_train-1)*8+8);
% XICs_C = XICs_C(:, (TrainingId_train-1)*8+1:(TrainingId_train-1)*8+8);
% IntervalListBig = IntervalListBig(TrainingId_train, :);
% INTERVAL_Matrix = INTERVAL_Matrix(TrainingId_train, :);

% TrainingId_AB = TrainingId_train(1:500);
% TrainingId_AC = TrainingId_train(1:500);
% TrainingId_BC = TrainingId_train(1:500);

% TrainingId_train = genTrainingId(INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,5)/2);
% TrainingId_AB = TrainingId_train(1:200);
% TrainingId_AC = genTrainingId(INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,9)/3);
% TrainingId_BC = genTrainingId(INTERVAL_Matrix(:,5)/2, INTERVAL_Matrix(:,9)/3);
% Diagnostic Code End------------------------------------------

% TrainingId_AB = INTERVAL_Matrix_TrainingSet(INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,5)/2);
% TrainingId_AC = genTrainingId(INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,9)/3);
% TrainingId_BC = genTrainingId(INTERVAL_Matrix(:,5)/2, INTERVAL_Matrix(:,9)/3);


TrainingId_ABTmp = genTrainingId(INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,5)/2);
TrainingId_ACTmp = genTrainingId(INTERVAL_Matrix(:,1)/1, INTERVAL_Matrix(:,9)/3);
TrainingId_BCTmp = genTrainingId(INTERVAL_Matrix(:,5)/2, INTERVAL_Matrix(:,9)/3);

TrainingId_AB = TrainingId_ABTmp(1:120);
TrainingId_AC = TrainingId_ACTmp(1:80);
TrainingId_BC = TrainingId_BCTmp(1:5);

% Diagnostic code starts
% load TrId
% TrainingId_AB = TrId;
% TrainingId_AC = TrId;
% TrainingId_BC = TrId;
% Diagnostic code ends

%  trainingModelIndicator = [length(TrainingId_AB)>100; length(TrainingId_AC)>100; length(TrainingId_BC)>100];
%%%%%%%%%%% build models AT AR AKL (AB AC BC datasets) and save them in
%%%%%%%%%%% Model_Training with AB AC BC order

%%%%%%%%%%%%%%%%%%%%%%%Generate Warping function of AB AC BC
PP_Parameter.PP_AB = genPP_parameter(retentiont_A, retentiont_B, Information_Matrix_Total, TrainingId_AB, 1);
PP_Parameter.PP_AC = genPP_parameter(retentiont_A, retentiont_C, Information_Matrix_Total, TrainingId_AC, 2);
PP_Parameter.PP_BC = genPP_parameter(retentiont_B, retentiont_C, Information_Matrix_Total, TrainingId_BC, 3);

%%%%%%%%%%%%%%%%%%%%%%%
[TrainingInformation_AB, Training_Information_Matrix_AB] = trainModel(retentiontl1_A, retentiont_A, XICs_A, retentiontl1_B, retentiont_B, XICs_B, Information_Matrix_Total, IntervalListBig, INTERVAL_Matrix, TrainingId_AB, 1);
[TrainingInformation_AC, Training_Information_Matrix_AC] = trainModel(retentiontl1_A, retentiont_A, XICs_A, retentiontl1_C, retentiont_C, XICs_C, Information_Matrix_Total, IntervalListBig, INTERVAL_Matrix, TrainingId_AC, 2);
[TrainingInformation_BC, Training_Information_Matrix_BC] = trainModel(retentiontl1_B, retentiont_B, XICs_B, retentiontl1_C, retentiont_C, XICs_C, Information_Matrix_Total, IntervalListBig, INTERVAL_Matrix, TrainingId_BC, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model_Training_AB = Generate_Training_Scores(TrainingInformation_AB, PP_Parameter, 1);
Model_Training_AC = Generate_Training_Scores(TrainingInformation_AC, PP_Parameter, 2);
Model_Training_BC = Generate_Training_Scores(TrainingInformation_BC, PP_Parameter, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters{1} =BuildATARAKLModelParameters(TrainingInformation_AB,Model_Training_AB);
Parameters{2} =BuildATARAKLModelParameters(TrainingInformation_AC,Model_Training_AC);
Parameters{3} =BuildATARAKLModelParameters(TrainingInformation_BC,Model_Training_BC);
%%%%%%%%%%%%%%%%%%%%%%%%