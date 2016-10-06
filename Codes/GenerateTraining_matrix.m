function Training_matrix=GenerateTraining_matrix(Training_candidate_matrixv1,...
    groundtruthinterval01_final,groundtruthinterval02_final,groundtruthinterval03_final,...
    IntervalList01_final,IntervalList02_final,IntervalList03_final,...
    monoXICs01_final,monoXICs02_final,monoXICs03_final)
Training_matrix=[];
for i=1:size(Training_candidate_matrixv1,1)
    
    id_in01=Training_candidate_matrixv1(i,8);
    id_in02=Training_candidate_matrixv1(i,16);
    id_in03=Training_candidate_matrixv1(i,24);
    
    id_intervalindex01=groundtruthinterval01_final(id_in01,2);
    id_intervalindex02=groundtruthinterval02_final(id_in02,2);
    id_intervalindex03=groundtruthinterval03_final(id_in03,2);
    
    interval01=IntervalList01_final{id_in01}.intervallist(id_intervalindex01,:);
    interval02=IntervalList02_final{id_in02}.intervallist(id_intervalindex02,:);
    interval03=IntervalList03_final{id_in03}.intervallist(id_intervalindex03,:);
    
    Ground_elution_peak01=monoXICs01_final(interval01(1):interval01(2),id_in01);
    Ground_elution_peak02=monoXICs02_final(interval02(1):interval02(2),id_in02);
    Ground_elution_peak03=monoXICs03_final(interval03(1):interval03(2),id_in03);
    
    highest_peak01=max(Ground_elution_peak01);
    highest_peak02=max(Ground_elution_peak02);
    highest_peak03=max(Ground_elution_peak03);
    
    if  highest_peak01>=1000000 && highest_peak02>=1000000 && highest_peak03>=1000000
        Training_matrix=[Training_matrix;Training_candidate_matrixv1(i,:)];        
    end
end
