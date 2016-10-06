function [intervalsdata01v2 ,Scan_start_afterfilter,Scan_end_afterfilter]=MS2intervalfilter_SILAC(intervaldata, MS2_Peak_Intensity01, intervalsdata01v1_afterfilter, MS2_Peak_Intensity01_afterfilter, Scan_start, Scan_end, Ms2_scannumber, threshold)
    
   
    Spearman_Corr01=zeros(1,size(intervaldata,1));
    Spearman_Corr01_afterfilter=zeros(1,size(intervaldata,1));
    Top_isotope_num=8; %%%% O18 6 isotope position ; SILAC 8 isotope position
%     plot_matirx(intervaldata,6);
%     plot_matirx(intervaldata,8);
%     plot_matirx(intervalsdata01v1_afterfilter,6);
%     plot_matirx(intervalsdata01v1_afterfilter,8); 
    
    for k=1:size(intervaldata,1)
        [V_ms2_Int01,P_ms2_Int01]=sort(MS2_Peak_Intensity01,'descend');
        a_toppeak=intervaldata(k,P_ms2_Int01(1:Top_isotope_num));
        Spearman_Corr01_toppeak(k)=corr(a_toppeak',MS2_Peak_Intensity01(P_ms2_Int01(1:Top_isotope_num))','type','spearman');
    
        a=intervaldata(k,:);
        Spearman_Corr01(k)=corr(a',MS2_Peak_Intensity01','type','spearman');    
        
        a_afterfilter=intervalsdata01v1_afterfilter(k,:);
        Spearman_Corr01_afterfilter(k)=corr(a_afterfilter',MS2_Peak_Intensity01_afterfilter','type','spearman'); 
        
        InterPeak_KL_distance(k)=KL_calculate(a_afterfilter,MS2_Peak_Intensity01_afterfilter); 

    end
    Goodpoint_p01=InterPeak_KL_distance<=threshold;
%     Goodpoint_p01=InterPeak_KL_distance<=-2.5;
%     Goodpoint_p01=Spearman_Corr01_afterfilter>=0.85;
%     Goodpoint_p01=Spearman_Corr01>=0.5;
%     Goodpoint_p01=Spearman_Corr01_toppeak>=0.5;
    
    Goodpoint_p01=[0,Goodpoint_p01,0];
    Rownumber01=Ms2_scannumber-Scan_start+1+1;
    Rownumber01_start=Rownumber01;
    while Goodpoint_p01(Rownumber01_start)==1
        Rownumber01_start=Rownumber01_start-1;
    end
    Rownumber01_end=Rownumber01;
    while Goodpoint_p01(Rownumber01_end)==1
        Rownumber01_end=Rownumber01_end+1;
    end
    Goodinterval_start01=Rownumber01_start;
    Goodinterval_end01=Rownumber01_end-2;
    intervalsdata01v2=intervaldata(Goodinterval_start01:Goodinterval_end01,:);
    Scan_start_afterfilter=Scan_start+Goodinterval_start01-1;
    Scan_end_afterfilter=Scan_start+Goodinterval_end01-1;