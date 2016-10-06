function [Parameters,Threshold]=BuildATModelParameters(Training_Information_Matrix,Model_Training)


 %%%%%%%%%%%%%%%%%% file count 1,2 and 3 is for model AB AC BC
 %%%%%%%%%%%%%%%%%%

    Corresponding_AT=[];Corresponding_AR=[];
    NonCorresponding_AT=[];NonCorresponding_AR=[];
    
    for i=1:length(Model_Training)
        

        MS2_ID01=Training_Information_Matrix{i}{1}.MS2_ID;
        MS2_ID02=Training_Information_Matrix{i}{2}.MS2_ID;
            

        AT_Matrix=Model_Training{i}.AT;

        Corresponding_AT(i)=AT_Matrix(1,MS2_ID02);
        
        %V1=AT_Matrix(:,MS2_ID02);V1(MS2_ID01)=[];
        V2=AT_Matrix';V2(MS2_ID02)=[];
        NonCorresponding_AT=[NonCorresponding_AT; V2];
     
           
    end
    
    [ATmodel_MU_sig,ATmodel_SIGMA_sig] = normfit(Corresponding_AT);
    [ATmodel_MU_notsig,ATmodel_SIGMA_notsig] = normfit(NonCorresponding_AT);


    %%%%%%%%%%%%%% AT and log2(AKL) are guassian models with mu and sigma 1-2 for
    %%%%%%%%%%%%%% corresponding and 3-4 for non-corresponding
    %%%%%%%%%%%%%% AR is gamma model with two parameters organized as a
    %%%%%%%%%%%%%% 1by2 vector
    
Parameters.Para_ATmodel=[ATmodel_MU_sig,ATmodel_SIGMA_sig,ATmodel_MU_notsig,ATmodel_SIGMA_notsig];
    
PAT_corr_param = pdf('Normal',Corresponding_AT,ATmodel_MU_sig,ATmodel_SIGMA_sig);
PAT_noncorr_param = pdf('Normal',Corresponding_AT,ATmodel_MU_notsig,ATmodel_SIGMA_notsig);
Index_PATn_PATc=PAT_noncorr_param./PAT_corr_param;L_Before=length(Index_PATn_PATc);

% PAT_corr_param1 = pdf('Normal',NonCorresponding_AT,ATmodel_MU_sig,ATmodel_SIGMA_sig);
% PAT_noncorr_param1 = pdf('Normal',NonCorresponding_AT,ATmodel_MU_notsig,ATmodel_SIGMA_notsig);
% Index_PATn_PATc1=PAT_noncorr_param1./PAT_corr_param1;

InfPos=((Index_PATn_PATc)<=5);Index_PATn_PATc=Index_PATn_PATc(InfPos);L_After=length(Index_PATn_PATc);
Pecent_AT_Treshold=L_After/L_Before*100;


%Mean_PAT=mean(Index_PATn_PATc);VAR_PAT=var(Index_PATn_PATc);
Threshold=max(Index_PATn_PATc);%Mean_PAT+30*VAR_PAT;