function         [N_time_diff,LE_diff,TP_value_vector,LEP_value_vector]=Calculate_PV_TLEnulldistribution(TimeInterval_A,Interval_data_A,...
            isoList,TOF_retentiontl1,Orbit_Pep_N_time,Pep_labeling_eff,...
            PP_null,T_null_PARMHAT,R_null_PARMHAT,LE_null_PARMHAT)
        
        minf=0.0;maxf=1;rangerate=0.6;
        for i=1:length(Interval_data_A)
            
%             Elu_Profile_A=sum(Interval_data_A{i},2)/6;
%             Elu_Profile_B=sum(Interval_data_B{i},2)/6;
%             
%             [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(Elu_Profile_A', Elu_Profile_B');
%             [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
%             TOF_Rstatistic(k,1)=STATS(1); 
            
            NTime_A=sum(TimeInterval_A(i,:))/2/max(TOF_retentiontl1);
            N_time_diff(i,1)=abs(polyval(PP_null,NTime_A)-Orbit_Pep_N_time);
            
            intensity_A=sum(Interval_data_A{i},1);
            est=O18rateLinear(isoList(1:6),intensity_A,minf,maxf,rangerate,1);
            TOF_Pep_A_labeling_eff=est(3);
            TOF_Pep_A_labeling_O18rate=est(2);
            LE_diff(i,1)=abs(TOF_Pep_A_labeling_eff-Pep_labeling_eff);
            
            TP_value_vector(i)=gamcdf(N_time_diff(i,1),T_null_PARMHAT(1),T_null_PARMHAT(2));
            LEP_value_vector(i)=gamcdf(LE_diff(i,1),LE_null_PARMHAT(1),LE_null_PARMHAT(2));
            
        end
        
        
        
        
        