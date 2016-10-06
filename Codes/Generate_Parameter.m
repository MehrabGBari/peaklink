function [AR0102_PARMHAT1,AR0102_PARMHAT2,AR0201_PARMHAT1,AR0201_PARMHAT2,...
AR0103_PARMHAT1,AR0103_PARMHAT2,AR0301_PARMHAT1,AR0301_PARMHAT2,...
AR0203_PARMHAT1,AR0203_PARMHAT2,AR0302_PARMHAT1,AR0302_PARMHAT2,...
timemodel0102_MU_sig,timemodel0102_SIGMA_sig,...
timemodel0201_MU_sig,timemodel0201_SIGMA_sig,...
timemodel0103_MU_sig,timemodel0103_SIGMA_sig,...
timemodel0301_MU_sig,timemodel0301_SIGMA_sig,...
timemodel0203_MU_sig,timemodel0203_SIGMA_sig,...
timemodel0302_MU_sig,timemodel0302_SIGMA_sig,...
timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig,...
timemodel0201_MU_notsig,timemodel0201_SIGMA_notsig,...
timemodel0103_MU_notsig,timemodel0103_SIGMA_notsig,...
timemodel0301_MU_notsig,timemodel0301_SIGMA_notsig,...
timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig,...
timemodel0302_MU_notsig,timemodel0302_SIGMA_notsig]=Generate_Parameter(AT0102diff_corr,AT0102diff_non_corr,...
    AT0201diff_corr,AT0201diff_non_corr,...
    AT0103diff_corr,AT0103diff_non_corr,...
    AT0301diff_corr,AT0301diff_non_corr,...
    AT0203diff_corr,AT0203diff_non_corr,...
    AT0302diff_corr,AT0302diff_non_corr,...
    PR01diff_non_corr,PR01diff_corr,...
    PR02diff_non_corr,PR02diff_corr,...
    PR03diff_non_corr,PR03diff_corr,...
    AR0102diff_non_corr,AR0201diff_non_corr,AR0102diff_corr,...
    AR0103diff_non_corr,AR0301diff_non_corr,AR0103diff_corr,...
    AR0203diff_non_corr,AR0302diff_non_corr,AR0203diff_corr)

    [timemodel0102_MU_sig,timemodel0102_SIGMA_sig] = normfit(AT0102diff_corr);
    [timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig] = normfit(AT0102diff_non_corr);
    [timemodel0201_MU_sig,timemodel0201_SIGMA_sig] = normfit(AT0201diff_corr);
    [timemodel0201_MU_notsig,timemodel0201_SIGMA_notsig] = normfit(AT0201diff_non_corr);
    [timemodel0103_MU_sig,timemodel0103_SIGMA_sig] = normfit(AT0103diff_corr);
    [timemodel0103_MU_notsig,timemodel0103_SIGMA_notsig] = normfit(AT0103diff_non_corr);
    [timemodel0301_MU_sig,timemodel0301_SIGMA_sig] = normfit(AT0301diff_corr);
    [timemodel0301_MU_notsig,timemodel0301_SIGMA_notsig] = normfit(AT0301diff_non_corr);
    [timemodel0203_MU_sig,timemodel0203_SIGMA_sig] = normfit(AT0203diff_corr);
    [timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig] = normfit(AT0203diff_non_corr);
    [timemodel0302_MU_sig,timemodel0302_SIGMA_sig] = normfit(AT0302diff_corr);
    [timemodel0302_MU_notsig,timemodel0302_SIGMA_notsig] = normfit(AT0302diff_non_corr);

    aa= PR01diff_non_corr>=0 & PR01diff_non_corr<=1;
    PR01diff_non_corr_noNaN=PR01diff_non_corr(aa);
    aa= PR01diff_corr>=0 & PR01diff_corr<=1;
    PR01diff_corr_noNaN=PR01diff_corr(aa);
    PR01_PARMHAT1=gamfit(PR01diff_corr_noNaN);
    PR01_PARMHAT2=gamfit(PR01diff_non_corr_noNaN);

    aa= PR02diff_non_corr>=0 & PR02diff_non_corr<=1;
    PR02diff_non_corr_noNaN=PR02diff_non_corr(aa);
    aa= PR02diff_corr>=0 & PR02diff_corr<=1;
    PR02diff_corr_noNaN=PR02diff_corr(aa);
    PR02_PARMHAT1=gamfit(PR02diff_corr_noNaN);
    PR02_PARMHAT2=gamfit(PR02diff_non_corr_noNaN);
    
    aa= PR03diff_non_corr>=0 & PR03diff_non_corr<=1;
    PR03diff_non_corr_noNaN=PR03diff_non_corr(aa);
    aa= PR03diff_corr>=0 & PR03diff_corr<=1;
    PR03diff_corr_noNaN=PR03diff_corr(aa);
    PR03_PARMHAT1=gamfit(PR03diff_corr_noNaN);
    PR03_PARMHAT2=gamfit(PR03diff_non_corr_noNaN);

    aa= AR0102diff_corr>=0 & AR0102diff_corr<=1;
    AR0102diff_corr_noNaN=AR0102diff_corr(aa);
    aa= AR0102diff_non_corr>=0 & AR0102diff_non_corr<=1;
    AR0102diff_non_corr_noNaN=AR0102diff_non_corr(aa);
    aa= AR0201diff_non_corr>=0 & AR0201diff_non_corr<=1;
    AR0201diff_non_corr_noNaN=AR0201diff_non_corr(aa);
    AR0102_PARMHAT1=gamfit(1-AR0102diff_corr_noNaN);
    AR0102_PARMHAT2=gamfit(AR0102diff_non_corr_noNaN);
    AR0201_PARMHAT1=AR0102_PARMHAT1;
    AR0201_PARMHAT2=gamfit(AR0201diff_non_corr_noNaN);
    
    aa= AR0103diff_corr>=0 & AR0103diff_corr<=1;
    AR0103diff_corr_noNaN=AR0103diff_corr(aa);
    aa= AR0103diff_non_corr>=0 & AR0103diff_non_corr<=1;
    AR0103diff_non_corr_noNaN=AR0103diff_non_corr(aa);
    aa= AR0301diff_non_corr>=0 & AR0301diff_non_corr<=1;
    AR0301diff_non_corr_noNaN=AR0301diff_non_corr(aa);
    AR0103_PARMHAT1=gamfit(1-AR0103diff_corr_noNaN);
    AR0103_PARMHAT2=gamfit(AR0103diff_non_corr_noNaN);
    AR0301_PARMHAT1=AR0103_PARMHAT1;
    AR0301_PARMHAT2=gamfit(AR0301diff_non_corr_noNaN);
    
    aa= AR0203diff_corr>=0 & AR0203diff_corr<=1;
    AR0203diff_corr_noNaN=AR0203diff_corr(aa);
    aa= AR0203diff_non_corr>=0 & AR0203diff_non_corr<=1;
    AR0203diff_non_corr_noNaN=AR0203diff_non_corr(aa);
    aa= AR0302diff_non_corr>=0 & AR0302diff_non_corr<=1;
    AR0302diff_non_corr_noNaN=AR0302diff_non_corr(aa);
    AR0203_PARMHAT1=gamfit(1-AR0203diff_corr_noNaN);
    AR0203_PARMHAT2=gamfit(AR0203diff_non_corr_noNaN);
    AR0302_PARMHAT1=AR0203_PARMHAT1;
    AR0302_PARMHAT2=gamfit(AR0302diff_non_corr_noNaN);
    

