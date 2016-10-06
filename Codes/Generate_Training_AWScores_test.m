function [Result,AccuracyATAW1step_test]= Generate_Training_AWScores_test(TrainingInformation,PP_Parameter, flag,Parameters,Threshold)
%%%%%%%%%%%%%%%%%%% We need calculate AB AC BC scores of AT AR and AKL
%%%%%%%%%%%%%%%%%%% Totally need 3 loops
%%%%%%%%%%%%%%%%%%% Each contains two datasets a and b
switch flag
    
    case 1
        Index=[1,2];
        PP=PP_Parameter.PP_AB;
    case 2
        Index=[1,3];
        PP=PP_Parameter.PP_AC;
        
    case 3
        Index=[2,3];
        PP=PP_Parameter.PP_BC;
        
end

for i = 1:length(TrainingInformation)
    Training_a=TrainingInformation{i}{1};
    Training_b=TrainingInformation{i}{2};
    
    %%%%%%
    intervalmatrix01=Training_a.Interval_after7_combine;
    intervalmatrix02=Training_b.Interval_after7_combine;
    ID_ms2_a=Training_a.MS2_ID; % ms2 identified id
    ID_ms2_b=Training_b.MS2_ID; % ms2 identified id
    Elution_Profile_a=Training_a.IntervalElutionMatrix;
    Elution_Profile_b=Training_b.IntervalElutionMatrix;
    [v,p]=max(Training_a.iso);
    Elution_matrix=[Elution_Profile_a{ID_ms2_a};Elution_Profile_b{ID_ms2_b}];
    Vec=sum(Elution_matrix,1);
    if Vec(p)>=Vec(p+4)
        Id=p;
    else
        Id=p+4;
    end
    
    
    %for j1=1:size(intervalmatrix01,1)
 monoelutionprofile01=Elution_Profile_a{ID_ms2_a}(:,Id);

        for j2=1:size(intervalmatrix02,1)
             monoelutionprofile01=Elution_Profile_a{ID_ms2_a}(:,Id);

            monoelutionprofile02=Elution_Profile_b{j2}(:,Id);
            
               
            [V_a,P_a]=max(monoelutionprofile01);
            [V_b,P_b]=max(monoelutionprofile02);
            t1=Training_a.Interval_after7_combine_Time{ID_ms2_a}(P_a);
            t2=Training_b.Interval_after7_combine_Time{j2}(P_b);
            AT=polyval(PP,t1)-t2;
            PAT_cIndex = pdf('Normal',AT,Parameters.Para_ATmodel(1),Parameters.Para_ATmodel(2));
            PAT_ncIndex= pdf('Normal',AT,Parameters.Para_ATmodel(3),Parameters.Para_ATmodel(4));

            if (PAT_ncIndex/PAT_cIndex)<=Threshold
                
               monoelutionprofile01=Elution_Profile_a{ID_ms2_a};
               monoelutionprofile02=Elution_Profile_b{j2};

             for ii=1:8 
             LC_profile_A_CorrSum=Resampl(monoelutionprofile01(:,ii));
             LC_profile_B_CorrSum =Resampl(monoelutionprofile02(:,ii));
%              LC_profile_A_CorrSum=Resampl(monoelutionprofile01);
%              LC_profile_B_CorrSum =Resampl(monoelutionprofile02);
             [alignedProfile_A_Corr alignedProfile_B_Corr shortlength]=alignProfiles(LC_profile_A_CorrSum,LC_profile_B_CorrSum);
             [Corr_all Corr_9 Corr_7]=wavedec_corr(alignedProfile_A_Corr,alignedProfile_B_Corr);
             
            Corr_7Coeff(ii)=Corr_7;%3*Corr_7+Corr_9;%+Corr_all;
            
            
             end
             Corr_7=mean(Corr_7Coeff);

             
             AW(1,j2)= Corr_7/(PAT_ncIndex/PAT_cIndex) ;
             if Corr_7~=0
             PAT_c(1,j2)= PAT_cIndex; 
             PAT_nc(1,j2)= PAT_ncIndex;   
             else
             
             
             PAT_c(1,j2)= 0; 
             PAT_nc(1,j2)= 0; 
             end

             if j2==ID_ms2_b
                CorsCorr(1,i)=PAT_c(1,j2);
                CorsCorr(2,i)=PAT_nc(1,j2);
             end
            else
             
              AW(1,j2)=0;
              PAT_c(1,j2)=0; 
              PAT_nc(1,j2)=0; 
            end   
        end
 
        [V,P]=max(AW);
        if P==ID_ms2_b
           True(i)=1;
        else
          True(i)=0;   
        end
    %end
         Result{i}.AW=AW;
         Result{i}.PAT_c=PAT_c;
         Result{i}.PAT_nc=PAT_nc;
        clear  AW PAT_c PAT_nc Corr_7 CoRR
  
end
 AccuracyATAW1step_test=100*sum(True==1)/length(True)       
        
        
        
        

