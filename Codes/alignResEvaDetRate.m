function [score, nIntersect, nTesting, posFin] = alignResEvaDetRate(TrainingId_AB, TrainingId_AC, TestModel, Information_Matrix_Total, INTERVAL_Matrix)

nIntersect = zeros(1,2);
nTesting = zeros(2,3);
score = zeros(2, 3);
posFin = cell(1,2);
for flag = 1:2;
    switch flag
        case 1
            TrainingId = TrainingId_AB;
        case 2
            TrainingId = TrainingId_AC;
    end
    
    
    TestModel_forTest = TestModel(TrainingId);
    Information_Matrix_Total_forTest = Information_Matrix_Total(TrainingId,:);
    INTERVAL_Matrix_forTest = INTERVAL_Matrix(TrainingId,:);
    backGround = [INTERVAL_Matrix_forTest(:, 2), INTERVAL_Matrix_forTest(:, 6), INTERVAL_Matrix_forTest(:, 10)];
    
    
    dd = 1;
    lenTest = zeros(1, 3);
    for scoreType =[2,1,4];
        [FinalResult,  ~, ~, ~, ~, ~, pos, ~] = pickAlignedPairs(TestModel_forTest, Information_Matrix_Total_forTest, scoreType);
        
       
        
        pos_forTest = pos;
        FinalResult_forTest = FinalResult;
        
        
        testingID = find(pos_forTest>200);
        idTesting = zeros(length(FinalResult_forTest), 3);
        for i = 1:length(FinalResult_forTest)
            idTesting(i, :) = FinalResult_forTest{i}.Cal_Posi;
        end
        scoreTmp= sum(backGround(pos_forTest(testingID),:)== idTesting(testingID,:))/length(backGround(pos_forTest(testingID),:)); 
        score(flag, dd) = scoreTmp(flag+1);
        lenTest(dd) = length(testingID);
        dd = dd + 1;
    end
    posFin{flag} = pos;
    nIntersect(flag) = length(TrainingId);
    nTesting(flag,:) = lenTest;
end