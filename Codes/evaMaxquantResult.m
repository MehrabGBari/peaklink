function [hlrMean, proteinList, peptideList, rep1, rep2, rep3] = evaMaxquantResult(peptidepath1, peptidepath2, peptidepath3)

maxquant01 = readtext(peptidepath1,'\t');
maxquant02 = readtext(peptidepath2,'\t');
maxquant03 = readtext(peptidepath3,'\t');

peptideSeq01 = maxquant01(2:end, 1);
peptideSeq02 = maxquant02(2:end, 1);
peptideSeq03 = maxquant03(2:end, 1);

proteins01 = maxquant01(2:end, 28);
proteins02 = maxquant02(2:end, 28);
proteins03 = maxquant03(2:end, 28);

HLR01 = cell2mat(maxquant01(2:end, 35));
HLR02 = cell2mat(maxquant02(2:end, 35));
HLR03 = cell2mat(maxquant03(2:end, 35));

rep1.peptide = peptideSeq01;
rep2.peptide = peptideSeq02;
rep3.peptide = peptideSeq03;

rep1.protein = proteins01;
rep2.protein = proteins02;
rep3.protein = proteins03;

rep1.HLR = HLR01;
rep2.HLR = HLR02;
rep3.HLR = HLR03;

[tmp, ~, ~] = union(peptideSeq01, peptideSeq02);
[entryUnion, ~, ~] = union(tmp, peptideSeq03);


[~, Ia3, Ib_rep1] = intersect(entryUnion, peptideSeq01);
[~, Ia4, Ib_rep2] = intersect(entryUnion, peptideSeq02);
[~, Ia5, Ib_rep3] = intersect(entryUnion, peptideSeq03);

flagMat = zeros(length(entryUnion), 3);

flagMat(Ia3, 1) = 1;
flagMat(Ia4, 2) = 1;
flagMat(Ia5, 3) = 1;

hlrMat = zeros(length(entryUnion), 3);

hlrMat(Ia3, 1) = HLR01(Ib_rep1);
hlrMat(Ia4, 2) = HLR02(Ib_rep2);
hlrMat(Ia5, 3) = HLR03(Ib_rep3);

proteinMat = cell(length(entryUnion), 3);
proteinMat(Ia3, 1) = proteins01(Ib_rep1);
proteinMat(Ia4, 2) = proteins02(Ib_rep2);
proteinMat(Ia5, 3) = proteins03(Ib_rep3);

peptideMat = cell(length(entryUnion), 3);
peptideMat(Ia3, 1) = peptideSeq01(Ib_rep1);
peptideMat(Ia4, 2) = peptideSeq02(Ib_rep2);
peptideMat(Ia5, 3) = peptideSeq03(Ib_rep3);

hlrMean = zeros(length(entryUnion), 1);
proteinList = cell(length(entryUnion), 1);
peptideList = cell(length(entryUnion), 1);
for i = 1:length(entryUnion)
    pos = find(flagMat(i, :)==1);
    hlrMean(i) = mean(hlrMat(i, pos));
    proteinList{i} = proteinMat{i, pos(1)};
    peptideList{i} = peptideMat{i, pos(1)};
end



