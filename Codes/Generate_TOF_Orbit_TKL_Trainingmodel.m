function  [T_diff,TOF_KL_diff]=Generate_TOF_Orbit_TKL_Trainingmodel(Spec_pepseq,Interval_data,Orbit_NT,TOF_NT)
                  
    Pepseq=Spec_pepseq;
    [peptidenew, massdiffList, isHeavy]=modprocess({Pepseq});
    peptidenew=peptidenew{1};
    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(Pepseq));
    Theory_mono1stiso_pattern=isotopepattern(1:2)/sum(isotopepattern(1:2));

    for k=1:length(Interval_data)
        T_diff(k,:)=[TOF_NT(k),Orbit_NT,TOF_NT(k)-Orbit_NT];                        
        Sum_vector=sum(Interval_data{k},1);
        Peak_mono_1st_pattern_combine=Sum_vector(1:2)/sum(Sum_vector(1:2));
        %%%%%%%%%%%%% log KL
        TOF_KL_diff(k,1)=KL_calculate(Peak_mono_1st_pattern_combine,Theory_mono1stiso_pattern);
        %%%%%%%%%%%%%
    end                  