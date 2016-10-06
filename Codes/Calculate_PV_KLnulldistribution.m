function    [KL_value_vector,KLP_value_corr_vector,KLP_value_non_corr_vector,...
                    KLProb_corr_vector,KLProb_non_corr_vector]=Calculate_PV_KLnulldistribution(peptideseqence,Interval_data,...
                    KL_pdf_pp_corr,KL_pdf_pp_non_corr,KL_cdf_pp_corr,KL_cdf_pp_non_corr)
        
            Pepseq=peptideseqence;
            [peptidenew, massdiffList, isHeavy]=modprocess({Pepseq});
            peptidenew=peptidenew{1};
            [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(Pepseq));
            Theory_mono1stiso_pattern=isotopepattern(1:2)/sum(isotopepattern(1:2));

             for k=1:length(Interval_data)                    
                Sum_vector=sum(Interval_data{k},1);
                Peak_mono_1st_pattern_combine=Sum_vector(1:2)/sum(Sum_vector(1:2));
                %%%%%%%%%%%%% log KL
                KL_value_vector(k,1)=KL_calculate(Peak_mono_1st_pattern_combine,Theory_mono1stiso_pattern);
                %%%%%%%%%%%%%
                KLProb_corr_vector(k,1)=ppval(KL_pdf_pp_corr,KL_value_vector(k,1));
                KLProb_non_corr_vector(k,1)=ppval(KL_pdf_pp_non_corr,KL_value_vector(k,1));
                
                KLP_value_corr_vector(k,1)=ppval(KL_cdf_pp_corr,KL_value_vector(k,1));
                KLP_value_non_corr_vector(k,1)=ppval(KL_cdf_pp_non_corr,KL_value_vector(k,1));          
             end