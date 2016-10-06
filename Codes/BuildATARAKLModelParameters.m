function Parameters=BuildATARAKLModelParameters(Training_Information_Matrix,Model_Training)


 %%%%%%%%%%%%%%%%%% file count 1,2 and 3 is for model AB AC BC
 %%%%%%%%%%%%%%%%%%

    Corresponding_AT=[];Corresponding_AR=[];Corresponding_AKL=[];
    NonCorresponding_AT=[];NonCorresponding_AR=[];NonCorresponding_AKL=[];
    
    for i=1:length(Model_Training)
        

        MS2_ID01=Training_Information_Matrix{i}{1}.MS2_ID;
        MS2_ID02=Training_Information_Matrix{i}{2}.MS2_ID;
            

        AT_Matrix=Model_Training{i}.AT;
        AR_Matrix=Model_Training{i}.AR;
        AKL_Matrix=Model_Training{i}.AKL;

        Corresponding_AT(i)=AT_Matrix(MS2_ID01,MS2_ID02);
        Corresponding_AR(i)=AR_Matrix(MS2_ID01,MS2_ID02);
        Corresponding_AKL(i)=AKL_Matrix(MS2_ID01,MS2_ID02);
        
        V1=AT_Matrix(:,MS2_ID02);V1(MS2_ID01)=[];
        V2=AT_Matrix(MS2_ID01,:)';V2(MS2_ID02)=[];
        NonCorresponding_AT=[NonCorresponding_AT; V1; V2];
    
        V1=AR_Matrix(:,MS2_ID02);V1(MS2_ID01)=[];
        V2=AR_Matrix(MS2_ID01,:)';V2(MS2_ID02)=[];
        NonCorresponding_AR=[NonCorresponding_AR; V1; V2];        
        
        V1=AKL_Matrix(:,MS2_ID02);V1(MS2_ID01)=[];
        V2=AKL_Matrix(MS2_ID01,:)';V2(MS2_ID02)=[];
        NonCorresponding_AKL=[NonCorresponding_AKL; V1; V2];      
    end
    
    [ATmodel_MU_sig,ATmodel_SIGMA_sig] = normfit(Corresponding_AT);
    [ATmodel_MU_notsig,ATmodel_SIGMA_notsig] = normfit(NonCorresponding_AT);
    AR_PARMHAT_sig=gamfit(Corresponding_AR);
    AR_PARMHAT_notsig=gamfit(NonCorresponding_AR);
    [AKLmodel_MU_sig,AKLmodel_SIGMA_sig] = normfit(Corresponding_AKL);
    [AKLmodel_MU_notsig,AKLmodel_SIGMA_notsig] = normfit(NonCorresponding_AKL);

    %%%%%%%%%%%%%% AT and log2(AKL) are guassian models with mu and sigma 1-2 for
    %%%%%%%%%%%%%% corresponding and 3-4 for non-corresponding
    %%%%%%%%%%%%%% AR is gamma model with two parameters organized as a
    %%%%%%%%%%%%%% 1by2 vector
     
    Parameters.Para_ATmodel=[ATmodel_MU_sig,ATmodel_SIGMA_sig,ATmodel_MU_notsig,ATmodel_SIGMA_notsig];
    Parameters.Para_AKLmodel=[AKLmodel_MU_sig,AKLmodel_SIGMA_sig,AKLmodel_MU_notsig,AKLmodel_SIGMA_notsig];
    Parameters.Para_ARmodel=[AR_PARMHAT_sig,AR_PARMHAT_notsig];
    

