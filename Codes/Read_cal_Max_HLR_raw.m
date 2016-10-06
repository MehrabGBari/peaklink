function [Peptideseq,Resultmatrix]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio)
load([Result_Save_Path,'\Pepinformation'])
load([Result_Save_Path,'\IntervalList01']);
load([Result_Save_Path,'\Result_HLRatio_SumPo1']);

for i=1:size(posi01,1)
    Peptideseq{i}=pep01{posi01(i,1)};
end
Resultmatrix(:,1)=ones(length(Cal_HLRatio),1)*PredefindeRatio;
Resultmatrix(:,2)=PepCal_Maxquant_peptide_HLR';
Resultmatrix(:,3)=Cal_HLRatio';




