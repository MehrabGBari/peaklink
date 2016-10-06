function  [LE_diff,LEP_value_corr_vector,LEP_value_non_corr_vector,...
                    LEProb_corr_vector,LEProb_non_corr_vector]=Calculate_PV_LEnulldistribution(Interval_data_A,...
                    isoList,Pep_labeling_eff,...
                    LE_pdf_pp_corr,LE_pdf_pp_non_corr,LE_cdf_pp_corr,LE_cdf_pp_non_corr)
        
        minf=0.0;maxf=1;rangerate=0.6;
        for i=1:length(Interval_data_A) 
            
            intensity_A=sum(Interval_data_A{i},1);
            est=O18rateLinear(isoList(1:6),intensity_A,minf,maxf,rangerate,1);
            TOF_Pep_A_labeling_eff=est(3);
            TOF_Pep_A_labeling_O18rate=est(2);
            LE_diff(i,1)=TOF_Pep_A_labeling_eff-Pep_labeling_eff;
            
            LEProb_corr_vector(i,1)=ppval(LE_pdf_pp_corr,LE_diff(i,1));
            LEProb_non_corr_vector(i,1)=ppval(LE_pdf_pp_non_corr,LE_diff(i,1));
                
            LEP_value_corr_vector(i,1)=abs(ppval(LE_cdf_pp_corr,LE_diff(i,1))-ppval(LE_cdf_pp_corr,-LE_diff(i,1)));
            LEP_value_non_corr_vector(i,1)=abs(ppval(LE_cdf_pp_non_corr,LE_diff(i,1))-ppval(LE_cdf_pp_non_corr,-LE_diff(i,1)));    
            
        end
        
        
        
        
        