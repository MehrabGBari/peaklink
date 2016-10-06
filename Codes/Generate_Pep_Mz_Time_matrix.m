function [Int_Matrix_mz_time,mz_unique,Pep_mz,Interval_XIC]=Generate_Pep_Mz_Time_matrix(peakl,MZ_interval,Interval_start,Interval_end)

    Int_Matrix_mz_time_List=[];
    for i=Interval_start:Interval_end
        ID_inmzInterval=find(peakl.scan(i).peak(:,1)>=MZ_interval(1) & peakl.scan(i).peak(:,1)<=MZ_interval(2));
        Int_Matrix_mz_time_List=[Int_Matrix_mz_time_List; peakl.scan(i).peak(ID_inmzInterval,:) i*ones(length(ID_inmzInterval),1)];
    end
    
    if ~isempty(Int_Matrix_mz_time_List)
   
            mz_unique=unique(Int_Matrix_mz_time_List(:,1));
            Int_Matrix_mz_time=zeros(length(mz_unique),Interval_end-Interval_start+1);
            for i=1:length(mz_unique)
                I=find(Int_Matrix_mz_time_List(:,1)==mz_unique(i));
                for k=1:length(I)
                    Int=Int_Matrix_mz_time_List(I(k),2);
                    Scan_id=Int_Matrix_mz_time_List(I(k),3);            
                    Int_Matrix_mz_time(i,Scan_id-Interval_start+1)=Int;            
                end
            end
            LC_profile=sum(Int_Matrix_mz_time,1);
            Weight_LC=LC_profile./sum(LC_profile);
            Int_raw_Weight=Int_Matrix_mz_time*Weight_LC';
        %     Int_raw_Weight=sum(Int_Matrix_mz_time,2);
            Int_Weight=Int_raw_Weight/sum(Int_raw_Weight);
            Pep_mz=mz_unique'*Int_Weight;


            MZ_profile=sum(Int_Matrix_mz_time,2);
            [Y,I]=max(MZ_profile);
            Timeprofile_standard=Int_Matrix_mz_time(I,:);
            for kkk=1:size(Int_Matrix_mz_time,1)
                    monoelutionprofile01=Timeprofile_standard;
                    monoelutionprofile02=Int_Matrix_mz_time(kkk,:);
                    if sum(monoelutionprofile02)~=0
                        [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01, monoelutionprofile02);
                        [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                        R_time_mz_matrix(kkk)=STATS(1);
                    else
                        R_time_mz_matrix(kkk)=0;                
                    end
            end
            ID_back=find(R_time_mz_matrix(I:end)<0.7);
            if ~isempty(ID_back)
                MZ_end=I+ID_back(1)-2;
            else 
                MZ_end=length(R_time_mz_matrix);
            end
            ID_front=find(R_time_mz_matrix(1:I)<0.7);
            if ~isempty(ID_front)
                MZ_start=ID_front(end)+1;
            else 
                MZ_start=1;
            end
            Interval_XIC=[sum(Int_Matrix_mz_time(MZ_start:MZ_end,:),1)]';
    
    else
        Int_Matrix_mz_time=zeros(1,Interval_end-Interval_start+1);
        mz_unique=sum(MZ_interval)/2;
        Pep_mz=sum(MZ_interval)/2;
        Interval_XIC=[sum(Int_Matrix_mz_time(1,:),1)]';
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    