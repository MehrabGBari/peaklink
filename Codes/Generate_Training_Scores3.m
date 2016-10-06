function Result=Generate_Training_Scores(TrainingInformation,PP_Parameter, flag)
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
    
    for j1=1:size(intervalmatrix01,1)
           monoelutionprofile01=Elution_Profile_a{j1};
        for j2=1:size(intervalmatrix02,1)
%             monoelutionprofile01=Elution_Profile_a{j1}(:,Id);
%             monoelutionprofile02=Elution_Profile_b{j2}(:,Id);
             %monoelutionprofile02=Elution_Profile_b{j2};      
             %if sum(monoelutionprofile01)~=0 && sum(monoelutionprofile02)~=0
%                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01', monoelutionprofile02');
%                [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
%                AR(j1,j2)=STATS(1);
% % % % % % %  LC_profile_A_CorrSum=Resampl(monoelutionprofile01);
% % % % % % %   LC_profile_B_CorrSum =Resampl(monoelutionprofile02);
% % % % % % %             [alignedProfile_A_Corr alignedProfile_B_Corr shortlength]=alignProfiles(LC_profile_A_CorrSum,LC_profile_B_CorrSum);
% % % % % % %             [Corr_all Corr_9 Corr_7]=wavedec_corr(alignedProfile_A_Corr,alignedProfile_B_Corr);
% % % % % % %               AR(j1,j2)=Corr_7;
          monoelutionprofile02=Elution_Profile_b{j2};
            
    for ii=1:8 
            LC_profile_A_CorrSum=Resampl(monoelutionprofile01(:,ii));
            LC_profile_B_CorrSum =Resampl(monoelutionprofile02(:,ii));
           [alignedProfile_A_Corr alignedProfile_B_Corr shortlength]=alignProfiles(LC_profile_A_CorrSum,LC_profile_B_CorrSum);
            %[alignedProfile_A_Corr alignedProfile_B_Corr]=alignProfiles_ver2(LC_profile_A_CorrSum,LC_profile_B_CorrSum);
            [Corr_all Corr_9 Corr_7]=wavedec_corr(alignedProfile_A_Corr,alignedProfile_B_Corr);
     Corr_7Coeff(ii)=Corr_7;
     
            [V_a,P_a]=max(monoelutionprofile01(:,ii));
            [V_b,P_b]=max(monoelutionprofile02(:,ii));
            t1=Training_a.Interval_after7_combine_Time{j1}(P_a);
            t2=Training_b.Interval_after7_combine_Time{j2}(P_b);
            AT1(ii)=polyval(PP,t1)-t2;
            
    end
   AR(j1,j2)=mean(Corr_7Coeff);
    AT(j1,j2)=mean(AT1);
% % % % % % %                 
% % % % % % %               
%            else
%                AR(j1,j2)=0;
           
            
%             [V_a,P_a]=max(Elution_Profile_a{j1}(:,Id));
%             [V_b,P_b]=max(Elution_Profile_b{j2}(:,Id));
%             t1=Training_a.Interval_after7_combine_Time{j1}(P_a);
%             t2=Training_b.Interval_after7_combine_Time{j2}(P_b);
%             AT(j1,j2)=polyval(PP,t1)-t2;
            
            Vector01=sum(Elution_Profile_a{j1},1);
            Vector02=sum(Elution_Profile_b{j2},1);
            AKL(j1,j2)=KL_calculate(Vector01,Vector02);
            
        end
    end
    Result{i}.AT=AT;
    Result{i}.AR=AR;
    Result{i}.AKL=AKL;
    
    clear AT AR AKL
end

        
        
        
        
        

