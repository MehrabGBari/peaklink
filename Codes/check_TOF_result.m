clc
clear all


load D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile_5574_f4\TOF_align_Result

load D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile_5574_f4\TOF_Final_Result_AT

for i=1:length(TOF_align_Result)
    Align_pep{i}=TOF_align_Result{i}.peptideseqence(3:end-2);
end

for i=1:length(TOF_Final_Result)
    Iden_pep{i}=TOF_Final_Result{i}.peptideseqence(3:end-2);
end

[Pep_cm,IA,IB]=intersect(Align_pep,Iden_pep);
length(Pep_cm)
length(TOF_Final_Result)

Iden_pep{IB}
ID=[];
for i=1:length(IA)
    if TOF_align_Result{IA(i)}.Judge_Align==0
        ID=[ID;i];
    end
end
TOF_align_Result{IA(889)}
TOF_Final_Result{IB(889)}



