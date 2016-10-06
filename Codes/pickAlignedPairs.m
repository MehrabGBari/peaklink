function [FinalResult, peptide, protein, TestModel, Information_Matrix_Total, iso, posEXT, flagA] = pickAlignedPairs(TestModel, Information_Matrix_Total, ChooseScoreMatrix)

n = length(TestModel);
FinalResult = cell(1, n);
peptide = cell(n, 1);
protein = cell(n, 1);
flagA = zeros(n, 2);
iso = cell(n, 1);

dd = 1;
for i = 1:length(TestModel)
    MS2_IDMat = [TestModel{i}.ImportantInf{1}.MS2_ID, TestModel{i}.ImportantInf{2}.MS2_ID, TestModel{i}.ImportantInf{3}.MS2_ID]~=0;
    if (MS2_IDMat(1)>0&&MS2_IDMat(2)>0)||(MS2_IDMat(1)>0&&MS2_IDMat(3)>0)
        posEXT(dd) = i;
        dd = dd + 1;
    end
end

TestModel = TestModel(posEXT);
Information_Matrix_Total = Information_Matrix_Total(posEXT, :);


infoString{1} = 'Choose intervals with best AT scores!';
infoString{2} = 'Choose intervals with best AR scores!';
infoString{3} = 'Choose intervals with best AKL scores!';
infoString{4} = 'Choose intervals with best AT+AR scores!';
infoString{5} = 'Choose intervals with best AT+AKL scores!';
infoString{6} = 'Choose intervals with best AR+AKL scores!';
infoString{7} = 'Choose intervals with best AT+AR+AKL scores!';

disp(infoString{ChooseScoreMatrix})



n = length(TestModel);
for i=1:n
        FinalResult{i}.ImportantInf=TestModel{i}.ImportantInf;
        %%%%%%%%%%%%%%% ScoreMatrix 7 cell for AT AR AKL AT+AR AT+AKL
        %%%%%%%%%%%%%%% AR+AKL AT+AR+AKL
        %%%%%%%%%%%%%%% 01 for 1st data pair 02 for 2nd data pair
        
        %%%%% 1:AT 2:AR 3:AKL 4:AT+AR 5:AT+AKL
        %%%%%%%%%%%%%%% 6:AR+AKL 7:AT+AR+AKL
        J_Vec=TestModel{i}.J_Vec;
        J_Val=sum(J_Vec*[1 2 4]');
        switch J_Val
            case 0%%[0 0 0]
              FinalResult{i}.J_Vec=J_Vec;
              peptide{i} = '-';
              protein{i} = '-';
              iso{i} = 0;
            case 1 %%[1 0 0]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});  
                MS2_ID=TestModel{i}.ImportantInf{J_Vec}.MS2_ID;
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_ID};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(MS2_ID,:));
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi01};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(MS2_ID,:));
                
                if Max_val01>0
                    flagA(i, 1) = 1;
                end
                if Max_val02>0
                    flagA(i, 2) = 1;
                end
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi02};
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[MS2_ID,Max_posi01,Max_posi02];
                
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};
            case 2 %%[0 1 0]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});  
                MS2_ID=TestModel{i}.ImportantInf{J_Vec}.MS2_ID;
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_ID};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(:,MS2_ID));  
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi01};
                
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(MS2_ID,:));
                
                if Max_val01>0
                    flagA(i, 1) = 1;
                end
                if Max_val02>0
                    flagA(i, 1) = 1;
                end
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi02};
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[Max_posi01,MS2_ID,Max_posi02];
                peptide{i} = Information_Matrix_Total{i, 9};
                protein{i} = Information_Matrix_Total{i, 14};
                iso{i} = Information_Matrix_Total{i, 16};
            case 4 %%[0 0 1]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID=TestModel{i}.ImportantInf{J_Vec}.MS2_ID;
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_ID};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(:,MS2_ID));
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi01};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(:,MS2_ID));
                
                if Max_val01>0
                    flagA(i, 1) = 1;
                end
                if Max_val02>0
                    flagA(i, 1) = 1;
                end
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi02};
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[Max_posi01,Max_posi02,MS2_ID];
                peptide{i} = Information_Matrix_Total{i, 17};
                protein{i} = Information_Matrix_Total{i, 22};
                iso{i} = Information_Matrix_Total{i, 24};
            case 3 %%[1 1 0]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID01=TestModel{i}.ImportantInf{1}.MS2_ID;
                MS2_ID02=TestModel{i}.ImportantInf{2}.MS2_ID;
                
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_ID01};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(MS2_ID01,:));
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_ID02};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(MS2_ID02,:));
                
                if Max_val01>0
                    flagA(i, 1) = 1;
                end
                if Max_val02>0
                    flagA(i, 1) = 1;
                end
                
                if Max_val01>=Max_val02                
                    FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi01};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,MS2_ID02,Max_posi01];
                else                
                    FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi02};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,MS2_ID02,Max_posi02];
                end
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};
                
            case 5 %%[1 0 1]
                
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID01=TestModel{i}.ImportantInf{1}.MS2_ID;
                MS2_ID02=TestModel{i}.ImportantInf{3}.MS2_ID;
                
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_ID01};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(MS2_ID01,:));
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_ID02};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(:,MS2_ID02));
                
                if Max_val01>0
                    flagA(i, 1) = 1;
                end
                if Max_val02>0
                    flagA(i, 1) = 1;
                end          
                if Max_val01>=Max_val02                
                    FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi01};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,Max_posi01,MS2_ID02];
                else                
                    FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi02};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,Max_posi02,MS2_ID02];
                end
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};
                
            case 6 %%[0 1 1]
                
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID01=TestModel{i}.ImportantInf{2}.MS2_ID;
                MS2_ID02=TestModel{i}.ImportantInf{3}.MS2_ID;
                
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_ID01};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(:,MS2_ID01));
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_ID02};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(:,MS2_ID02));
     
                if Max_val01>0
                    flagA(i, 1) = 1;
                end
                if Max_val02>0
                    flagA(i, 1) = 1;
                end         
                
                if Max_val01>=Max_val02                
                    FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi01};
                    FinalResult{i}.Cal_Posi=[Max_posi01,MS2_ID01,MS2_ID02];
                else                
                    FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi02};
                    FinalResult{i}.Cal_Posi=[Max_posi02,MS2_ID01,MS2_ID02];
                end
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                
                peptide{i} = Information_Matrix_Total{i, 9};
                protein{i} = Information_Matrix_Total{i, 14};
                iso{i} = Information_Matrix_Total{i, 16};
                
            case 7 %%[1 1 1]
                flagA(i, :) = [1,1];
                MS2_IDA=TestModel{i}.ImportantInf{1}.MS2_ID;
                MS2_IDB=TestModel{i}.ImportantInf{2}.MS2_ID;
                MS2_IDC=TestModel{i}.ImportantInf{3}.MS2_ID;                
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_IDA};
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_IDB};
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_IDC};
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[MS2_IDA,MS2_IDB,MS2_IDC];      
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};     
        end
end