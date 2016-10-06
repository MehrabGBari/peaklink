function     [ARmatrix,AKLmatrix]=CalARAKLMatrix(RowElutionProf,ColElutionProf)
nR = length(RowElutionProf);
nC = length(ColElutionProf);

if sum(sum(RowElutionProf{1}))~=0 && sum(sum(ColElutionProf{1}))~=0
    for j1=1:length(RowElutionProf)
        for j2=1:length(ColElutionProf)
            KLValueRow=sum(RowElutionProf{j1},1);
            KLValueCol=sum(ColElutionProf{j2},1);
            AKLmatrix(j1,j2)=KL_calculate(KLValueRow,KLValueCol);
            
            [V,Id]=max(KLValueRow+KLValueCol);            
            monoelutionprofile01=RowElutionProf{j1}(:,Id);
            monoelutionprofile02=ColElutionProf{j2}(:,Id);
     % % %       if sum(monoelutionprofile01)~=0 && sum(monoelutionprofile02)~=0
% % %                 [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01', monoelutionprofile02');
% % %                 [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
% % %                 ARmatrix(j1,j2)=STATS(1);
% % %                 
% % % %                 ARmatrix(j1,j2) = anyfunctionyouwanttobuild( monoelutionprofile01', monoelutionprofile02')
% % % 
                       LC_profile_A_CorrSum=Resampl(monoelutionprofile01);
             LC_profile_B_CorrSum =Resampl(monoelutionprofile02);
            [alignedProfile_A_Corr alignedProfile_B_Corr shortlength]=alignProfiles(LC_profile_A_CorrSum,LC_profile_B_CorrSum);
            [Corr_all Corr_9 Corr_7]=wavedec_corr(alignedProfile_A_Corr,alignedProfile_B_Corr);
             if Corr_7>0  
            ARmatrix(j1,j2)=Corr_7;
            else 
                ARmatrix(j1,j2)=0;                    
            end

        end
    end
else  
    AKLmatrix=zeros(nR, nC);
    ARmatrix=zeros(nR, nC);    
end