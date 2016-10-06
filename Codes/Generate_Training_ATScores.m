function Result=Generate_Training_ATScores(TrainingInformation,PP_Parameter, flag)
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
          % monoelutionprofile01=Elution_Profile_a{j1};
          monoelutionprofile01=Elution_Profile_a{ID_ms2_a}(:,Id);
        for j2=1:size(intervalmatrix02,1)
            %monoelutionprofile01=Elution_Profile_a{j1}(:,Id);
            monoelutionprofile02=Elution_Profile_b{j2}(:,Id);
             
            
            [V_a,P_a]=max(monoelutionprofile01);
            [V_b,P_b]=max(monoelutionprofile02);
            t1=Training_a.Interval_after7_combine_Time{ID_ms2_a}(P_a);
            t2=Training_b.Interval_after7_combine_Time{j2}(P_b);
            AT(1,j2)=polyval(PP,t1)-t2;
          
            
        end
    %end
    Result{i}.AT=AT;
   
    clear AT 
end
    

        
        
        
        
        

