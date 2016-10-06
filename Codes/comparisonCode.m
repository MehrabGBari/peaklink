% Evaluation Program

peptidepath1='C:\Users\Long\Desktop\SuperSilac\tumor02\combined\txt\peptides.txt';
peptidepath2='C:\Users\Long\Desktop\SuperSilac\tumor03\combined\txt\peptides.txt';
peptidepath3='C:\Users\Long\Desktop\SuperSilac\tumor04\combined\txt\peptides.txt';

[hlrMean, proteinList, peptideList, rep1, rep2, rep3] = evaMaxquantResult(peptidepath1, peptidepath2, peptidepath3);
save maxquant hlrMean proteinList peptideList

protein_rep1 = rep1.protein;
protein_rep2 = rep2.protein;
protein_rep3 = rep3.protein;

HLR_rep1 = rep1.HLR;
HLR_rep2 = rep2.HLR;
HLR_rep3 = rep3.HLR;

peptide_rep1 = rep1.peptide;
peptide_rep2 = rep2.peptide;
peptide_rep3 = rep3.peptide;

n = 200;
LRD_rep1 = zeros(n, 1);
LRD_rep2 = zeros(n, 1);
LRD_rep3 = zeros(n, 1);

for i = 1:200
    LRD_rep1(i) = calUniqueLRD(protein_rep1, HLR_rep1);
end
mean(LRD_rep1)
for i = 1:200
    LRD_rep2(i) = calUniqueLRD(protein_rep2, HLR_rep2);
end
mean(LRD_rep2)
for i = 1:200
    LRD_rep3(i) = calUniqueLRD(protein_rep3, HLR_rep3);
end
mean(LRD_rep3)




n = length(peptideSel);
labeledAminNum = zeros(n, 1);
for i = 1:n
    pos_K = strfind(peptideSel{i}, 'K');
    pos_R = strfind(peptideSel{i}, 'R');
    labeledAminNum(i) = length(pos_K) + length(pos_R);
end

repIndex =3;

pos = find(labeledAminNum <= 1);
proteinSel2 = proteinSel(pos);
HLR2 = HLR(pos, :);

index01 = findUniquePep(proteinList);
index02 = findUniquePep(proteinSel2);

uniProtMaxquant = proteinList(index01);
uniProtCui = proteinSel2(index02);

uniHLRMax = hlrMean(index01);
uniHLRCui = mean(HLR2(index02,repIndex), 2);

pepIntersect = intersect(uniProtMaxquant, uniProtCui);

posMaxquant = findSharedProt(uniProtMaxquant, pepIntersect);
posCui = findSharedProt(uniProtCui, pepIntersect);

uniProtSharedMaxq = uniProtMaxquant(posMaxquant);
uniHLRMaxShared = uniHLRMax(posMaxquant);

uniProtSharedCui = uniProtCui(posCui);
uniHLRCuiShared = uniHLRCui(posCui);

numofSample = 200;
LRD1 = zeros(numofSample, 1);
for i = 1:numofSample
    LRD1(i) = calLRD(uniProtSharedMaxq, uniHLRMaxShared);
end


LRD2 = zeros(numofSample, 1);

for i = 1:numofSample
    LRD2(i) = calLRD(uniProtSharedCui, uniHLRCuiShared);
end

mean(LRD1)
mean(LRD2)
% 
% mu = mean(keepLRD2 - LRD2)
% sigma = std(keepLRD2 - LRD2)
% 
%  P = normcdf(0,mu,sigma) 
% 
% 
% keepLRD2 = LRD2

i =2;
MS2_ID = FinalResultSel{i}.ImportantInf{1}.MS2_ID
plot(FinalResultSel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_ID})





