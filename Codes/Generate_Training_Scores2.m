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
    h=[8 7 6 6 5 5 4 4 3 3 2 2 1 1 1 1 2 2 3 3 4 4 5 5 6 6 7 ];h=h/sum(h);    
    for j1=1:size(intervalmatrix01,1)
        for j2=1:size(intervalmatrix02,1)
            
            monoelutionprofile01=Elution_Profile_a{j1}(:,Id);
            monoelutionprofile02=Elution_Profile_b{j2}(:,Id);
            if sum(monoelutionprofile01)~=0 && sum(monoelutionprofile02)~=0

                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01', monoelutionprofile02');
                newdata_sampling01=conv(h,newdata_sampling01);
                newdata_sampling02=conv(h,newdata_sampling02);
                [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                AR(j1,j2)=STATS(1);
                
                
%               AR(j1,j2) = anyfunctionyouwanttobuild( monoelutionprofile01', monoelutionprofile02')
                
            else
                AR(j1,j2)=0;
            end
            
            [V_a,P_a]=max(monoelutionprofile01);
            [V_b,P_b]=max(monoelutionprofile02);
            t1=Training_a.Interval_after7_combine_Time{j1}(P_a);
            t2=Training_b.Interval_after7_combine_Time{j2}(P_b);
            AT(j1,j2)=polyval(PP,t1)-t2;
            
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

        
        
        
        
        

