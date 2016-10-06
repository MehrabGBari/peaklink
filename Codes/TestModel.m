function [TestingInformation] = TestModel(retentiontl1_A, retentiont_A, XICs_A, retentiontl1_B, retentiont_B, XICs_B, Information_Matrix_Total, IntervalListBig, INTERVAL_Matrix, Testing_ID, flag)


n = length(Testing_ID);
switch flag
    case 1
        indexR = [1, 2];
    case 2
        indexR = [1, 3];
    case 3
        indexR = [2, 3];
end

IsoList=cell2mat(Information_Matrix_Total(Testing_ID,indexR(1)*8));
MS2scan_number_vector=[cell2mat(Information_Matrix_Total(Testing_ID,4)), cell2mat(Information_Matrix_Total(Testing_ID,12)), cell2mat(Information_Matrix_Total(Testing_ID,20))];
%%%%%%%%%%%%%%%% combine peptide information of A B C before building
%%%%%%%%%%%%%%%% models

IntervalSelA = IntervalListBig(Testing_ID,indexR(1));
IntervalSelB = IntervalListBig(Testing_ID,indexR(2));

MS2_IDA=INTERVAL_Matrix(Testing_ID,(indexR(1)-1)*4+2);
MS2_IDB=INTERVAL_Matrix(Testing_ID,(indexR(2)-1)*4+2);

MS2_SC_A=MS2scan_number_vector(:, indexR(1));
MS2_SC_B=MS2scan_number_vector(:, indexR(2));

Testing_Information_Matrix = cell(n, 2);
TestingInformation = cell(n, 1);

for i = 1:n
    intervalA = IntervalSelA{i}.intervallist_after7_combine;
    intervalB = IntervalSelB{i}.intervallist_after7_combine;
    
    for Interval_Row=1:size(intervalA,1)
        IntervalElutionMatrixA{Interval_Row}=XICs_A(intervalA(Interval_Row,1):intervalA(Interval_Row,2),(Testing_ID(i)-1)*8+1:Testing_ID(i)*8);
        Interval_after7_combine_TimeA{Interval_Row}=retentiontl1_A(intervalA(Interval_Row,1):intervalA(Interval_Row,2));
    end
    
    for Interval_Row=1:size(intervalB,1)
        IntervalElutionMatrixB{Interval_Row}=XICs_B(intervalB(Interval_Row,1):intervalB(Interval_Row,2),(Testing_ID(i)-1)*8+1:Testing_ID(i)*8);
        Interval_after7_combine_TimeB{Interval_Row}=retentiontl1_B(intervalB(Interval_Row,1):intervalB(Interval_Row,2));
    end
    
    
    Testing_Information_Matrix{i,1}.iso=IsoList(i, :);
    Testing_Information_Matrix{i,1}.Interval_after7_combine=intervalA;
    Testing_Information_Matrix{i,1}.MS2_ID=MS2_IDA(i);
    Testing_Information_Matrix{i,1}.Interval_after7_combine_Time=Interval_after7_combine_TimeA;
    Testing_Information_Matrix{i,1}.MS2timePoint=retentiont_A(MS2_SC_A(i));
    Testing_Information_Matrix{i,1}.IntervalElutionMatrix=IntervalElutionMatrixA;
    
    Testing_Information_Matrix{i,2}.iso=IsoList(i, :);
    Testing_Information_Matrix{i,2}.Interval_after7_combine=intervalB;
    Testing_Information_Matrix{i,2}.MS2_ID=MS2_IDB(i);
    Testing_Information_Matrix{i,2}.Interval_after7_combine_Time=Interval_after7_combine_TimeB;
    Testing_Information_Matrix{i,2}.MS2timePoint=retentiont_B(MS2_SC_B(i));
    Testing_Information_Matrix{i,2}.IntervalElutionMatrix=IntervalElutionMatrixB;
    
    TestingInformation{i} = Testing_Information_Matrix(i, :);
end