function [Cal_HLRatio_intersect_final_Sumpo1,Maxquant_HLR_K_Labeled_intersect_final_Sumpo1,...
    Cal_HLRatio_intersect_final_Sumpo2,Maxquant_HLR_K_Labeled_intersect_final_Sumpo2]=Read_cal_Max_HLR(Result_Save_Path)
load([Result_Save_Path,'\Result_HLRatio_SumPo1_final'])
Cal_HLRatio_intersect_final_Sumpo1=Cal_HLRatio_intersect_final;
Maxquant_HLR_K_Labeled_intersect_final_Sumpo1=Maxquant_HLR_K_Labeled_intersect_final;
clear Cal_HLRatio_intersect_final Maxquant_HLR_K_Labeled_intersect_final
load([Result_Save_Path,'\Result_HLRatio_SumPo2_final'])
Cal_HLRatio_intersect_final_Sumpo2=Cal_HLRatio_intersect_final;
Maxquant_HLR_K_Labeled_intersect_final_Sumpo2=Maxquant_HLR_K_Labeled_intersect_final;