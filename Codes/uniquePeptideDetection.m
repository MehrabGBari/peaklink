function [index1_r1, index2_r1] = uniquePeptideDetection(proteinSingle)

[uniProteinSeq,idA,idB] = unique(proteinSingle);
nn = 1;
for i=1:length(idA)
    idUni{i} = find(idB == i);
    n = length(idUni{i});
    if n > 1
        Id_abund = randsample(1:n, n, 'false');
        index1_r1(nn) = idUni{i}(Id_abund(1));
        index2_r1(nn) = idUni{i}(Id_abund(2));
        clear abund
        nn = nn + 1;
    end
end

if strcmp(proteinSingle(index1_r1),proteinSingle(index2_r1)) ==1
    disp('Unique peptide detection pass!')
else
    disp('Unique peptde detection error, please check!')
end
