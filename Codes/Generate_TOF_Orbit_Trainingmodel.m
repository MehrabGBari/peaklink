function        [TOF_Rstatistic,TOF_Pep_labeling_eff,T_diff,TOF_KL_diff]=Generate_TOF_Orbit_Trainingmodel(TOF_align_Result,Orbit_IntervalData,Orbit_NT,TOF_NT)
                    Orbit_Total_Elution_Profile=sum(Orbit_IntervalData,2)/6;
                    
                    
                    Pepseq=TOF_align_Result.peptideseqence;
                    [peptidenew, massdiffList, isHeavy]=modprocess({Pepseq});
                    peptidenew=peptidenew{1};
                    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(Pepseq));
                    Theory_mono1stiso_pattern=isotopepattern(1:2)/sum(isotopepattern(1:2));
                    
                    for k=1:length(TOF_align_Result.Interval_data)
                        T_diff(k,:)=[TOF_NT(k),Orbit_NT,TOF_NT(k)-Orbit_NT];                        
                        TOF_IntervalData=TOF_align_Result.Interval_data{k};
                        isoList=TOF_align_Result.iso;
                        minf=0.0;maxf=1;rangerate=0.6;
                        intensity=sum(TOF_IntervalData,1);
                        est=O18rateLinear(isoList(1:6),intensity,minf,maxf,rangerate,1);
                        TOF_Pep_labeling_eff(k,1)=est(3);
                        TOF_Pep_labeling_O18rate(k,1)=est(2);
                        
                        TOF_Total_Elution_Profile=sum(TOF_IntervalData,2)/6;
                        
                        [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(Orbit_Total_Elution_Profile', TOF_Total_Elution_Profile');
                        [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                        TOF_Rstatistic(k,1)=STATS(1);
                        
                        Com_matrix=[TOF_align_Result.Interval_data{k};TOF_align_Result.Interval_data_aligned{k}];
                        Peak_profile=sum(Com_matrix(:,1:2),2);
                        Weight=Peak_profile./sum(Peak_profile);
                        Peak_combine=sum(diag(Weight)*Com_matrix,1);
%                         
%                         Sum_vector=sum(TOF_align_Result.Interval_data{k},1);
%                         Sum_vector_aligned=sum(TOF_align_Result.Interval_data_aligned{k},1);
%                         Weight=sum(Sum_vector(1:2))/(sum(Sum_vector(1:2))+sum(Sum_vector_aligned(1:2)));
%                         Weight_aligned=sum(Sum_vector_aligned(1:2))/(sum(Sum_vector(1:2))+sum(Sum_vector_aligned(1:2)));
%                         Peak_combine=Weight*Sum_vector+Weight_aligned*Sum_vector_aligned;
                        Peak_mono_1st_pattern_combine=Peak_combine(1:2)/sum(Peak_combine(1:2));
                        %%%%%%%%%%%%% log KL
                        TOF_KL_diff(k,1)=KL_calculate(Peak_mono_1st_pattern_combine,Theory_mono1stiso_pattern);
                        %%%%%%%%%%%%%
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    