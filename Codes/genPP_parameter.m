function para = genPP_parameter(retentiont_A, retentiont_B, Information_Matrix_Total, TId, flag)

switch flag
    case 1
        indexR = [1, 2];
    case 2
        indexR = [1, 3];
    case 3
        indexR = [2, 3];
end
        
MS2scan_number_vector=[cell2mat(Information_Matrix_Total(TId,4)), cell2mat(Information_Matrix_Total(TId,12)), cell2mat(Information_Matrix_Total(TId,20))];
%%%%%%%%%%%%%%%% combine peptide information of A B C before building
%%%%%%%%%%%%%%%% models

MS2_SC_A=MS2scan_number_vector(:, indexR(1));
MS2_SC_B=MS2scan_number_vector(:, indexR(2));

MS2_TimeA =retentiont_A(MS2_SC_A);
MS2_TimeB =retentiont_B(MS2_SC_B);

para=polyfit(MS2_TimeA, MS2_TimeB, 4);




