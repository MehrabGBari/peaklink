
    function [LOGAKL,AR,AT,Normal_AT,T01,Normal_T01,T02,Normal_T02]=GenerateRTKL(intervals01,intervals02,XICs01,XICs02,Retentiontime01,Retentiontime02,iso_spec_pep)
    %%%%%% intervals01 and intervals02 are matrix;
    %%%%%% XICs01 and XICs02 are cell files that with same number of intervals.
    %%%%%% Each of the cell has matrixs (n*6) 6 colums are for different
    %%%%%% iso mz
    for k1=1:size(intervals01,1)
        for k2=1:size(intervals02,1)
    
            Scan_start01=intervals01(k1,1);
            Scan_end01=intervals01(k1,2);
            Scan_start02=intervals02(k2,1);
            Scan_end02=intervals02(k2,2);
            for MZ_posi=1:6
                LC_Interval_Matrix01(:,MZ_posi)=XICs01{k1}(:,MZ_posi);
                LC_Interval_Matrix02(:,MZ_posi)=XICs02{k2}(:,MZ_posi);
                if sum(LC_Interval_Matrix01(:,MZ_posi))~=0 && sum(LC_Interval_Matrix02(:,MZ_posi))~=0
                    monoelutionprofile01=LC_Interval_Matrix01(:,MZ_posi);
                    monoelutionprofile02=LC_Interval_Matrix02(:,MZ_posi);
                    [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01', monoelutionprofile02');
                    [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                    AR_iso_mz(MZ_posi)=STATS(1);
                else
                    AR_iso_mz(MZ_posi)=0;
                end
            end
            
            Vector01=sum(LC_Interval_Matrix01,1);
            Vector02=sum(LC_Interval_Matrix02,1);
            LOGAKL(k1,k2)=KL_calculate(Vector01,Vector02);

            if iso_spec_pep(1)>=iso_spec_pep(2)
                if AR_iso_mz(1)>=AR_iso_mz(5)
                    AR(k1,k2)=AR_iso_mz(1);
                else
                    AR(k1,k2)=AR_iso_mz(5);
                end
            else
                if AR_iso_mz(2)>=AR_iso_mz(6)
                    AR(k1,k2)=AR_iso_mz(2);
                else
                    AR(k1,k2)=AR_iso_mz(6);
                end
            end
            
            AT(k1,k2)=(Retentiontime01(Scan_start01)+Retentiontime01(Scan_end01))/2-(Retentiontime02(Scan_start02)+Retentiontime02(Scan_end02))/2;         
            Normal_AT(k1,k2)=(Retentiontime01(Scan_start01)+Retentiontime01(Scan_end01))/(2*max(Retentiontime01))-(Retentiontime02(Scan_start02)+Retentiontime02(Scan_end02))/(2*max(Retentiontime02));         
            T01(k1,k2)=(Retentiontime01(Scan_start01)+Retentiontime01(Scan_end01))/2;         
            Normal_T01(k1,k2)=(Retentiontime01(Scan_start01)+Retentiontime01(Scan_end01))/(2*max(Retentiontime01));         
            T02(k1,k2)=(Retentiontime02(Scan_start02)+Retentiontime02(Scan_end02))/2;         
            Normal_T02(k1,k2)=(Retentiontime02(Scan_start02)+Retentiontime02(Scan_end02))/(2*max(Retentiontime02));         
                          
            clear LC_Interval_Matrix01 LC_Interval_Matrix02
        end        
    end
    
    