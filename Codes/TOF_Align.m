function [TOF_align_Result,Com_ID_0102,IntervalList01,IntervalList02,...
    IntervalList_com01,IntervalList_com02]=TOF_Align(Pepcommon0102v1,TOF_mass,iso,TOF_charge,Protein_matrix_final,...
    TOF_XICs01,TOF_XICs02,...
    TOF_retentiont01l1,TOF_retentiont02l1)

pep01=Pepcommon0102v1;
posi01v1=[];
posi01v2=[];
posi01=[];
for i=1:length(pep01)
    TOF_XICs01_sepc_pep=zeros(size(TOF_XICs01{1},1),6);
    for j=1:6
        TOF_XICs01_sepc_pep(:,j)=TOF_XICs01{j}(:,i);        
    end
    [IntervalList,Jud_detect_good]=TOF_Verification_intervaldetection_Spec_Pep(pep01{i},iso(i,:),TOF_XICs01_sepc_pep);
    IntervalList01{i}=IntervalList;
    if Jud_detect_good==1
        posi01v1=[posi01v1;i];
        if IntervalList.intervallist_after7_combine(1,1)~=0 && IntervalList.intervallist_after7_combine(1,2)~=0
            posi01v2=[posi01v2;i];
% %             TOF_retentiont01l1v1=TOF_retentiont01l1';
%             T_m=TOF_retentiont01l1(IntervalList.intervallist_after7_combine)-IntervalList.ms2time;
%             T_v=T_m(:,1).*T_m(:,2);
%             Id_ms2interval=find(T_v<=0);
%             if ~isempty(Id_ms2interval)
%                 posi01=[posi01;i];
%             end            
        end
    end    
end

pep02=Pepcommon0102v1;
posi02v1=[];
posi02v2=[];
posi02=[];
for i=1:length(pep02)
    TOF_XICs02_sepc_pep=zeros(size(TOF_XICs02{1},1),6);
    for j=1:6
        TOF_XICs02_sepc_pep(:,j)=TOF_XICs02{j}(:,i);        
    end
    [IntervalList,Jud_detect_good]=TOF_Verification_intervaldetection_Spec_Pep(pep02{i},iso(i,:),TOF_ms2information02(i),TOF_XICs02_sepc_pep,TOF_retentiont02l1);
    IntervalList02{i}=IntervalList;
    if Jud_detect_good==1
        posi02v1=[posi02v1;i];
        if IntervalList.intervallist_after7_combine(1,1)~=0 && IntervalList.intervallist_after7_combine(1,2)~=0
            posi02v2=[posi02v2;i];
% %             TOF_retentiont02l1v1=TOF_retentiont02l1';
%             T_m=TOF_retentiont02l1(IntervalList.intervallist_after7_combine)-IntervalList.ms2time;
%             T_v=T_m(:,1).*T_m(:,2);
%             Id_ms2interval=find(T_v<=0);
%             if ~isempty(Id_ms2interval)
%                 posi02=[posi02;i];
%             end            
        end
    end    
end

save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\intervallist IntervalList01 IntervalList02
save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\IDs posi01v1 posi01v2 posi01 posi02v1 posi02v2 posi02

load D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\intervallist
load D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\IDs
% length(posi01v1)
% length(posi01v2)
% length(posi01)
% length(posi02v1)
% length(posi02v2)
% length(posi02)

%%%%%%%%%%%%%%%
[Com_ID_0102,ID_0102_01,ID_0102_02]=intersect(posi01v2,posi02v2);

for i=1:6
    TOF_XICs_com01{i}=TOF_XICs01{i}(:,posi01v2(ID_0102_01));
    TOF_XICs_com02{i}=TOF_XICs02{i}(:,posi02v2(ID_0102_02));
end
IntervalList_com01=IntervalList01(posi01v2(ID_0102_01));
IntervalList_com02=IntervalList02(posi02v2(ID_0102_02));

Pep_com=Pepcommon0102v1(Com_ID_0102);
% TOF_ms2information_com=TOF_ms2information01(Com_ID_0102);
iso_com=iso(Com_ID_0102,:);
TOF_mass_com=TOF_mass(Com_ID_0102);
TOF_charge_com=TOF_charge(Com_ID_0102);
Protein_matrix_final_com=Protein_matrix_final(Com_ID_0102);

save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\TOF_Total_information_cent0 Pep_com iso_com TOF_XICs_com01 TOF_XICs_com02 IntervalList_com01 IntervalList_com02   
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% check interval detection
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
Times_noise_std=3;%6;
instrument_h_int=10^7;
de_range=10^5;%10^4;
colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
for i=36:40%%length(posi01v2)
    ID_into=Com_ID_0102(i);
    figure
    subplot(2,1,1)
    for j=1:6
        th(j)=TOF_getThreshold_range(TOF_XICs01{j}(:,ID_into),instrument_h_int, de_range,Times_noise_std);
        plot(TOF_retentiont01l1,TOF_XICs01{j}(:,ID_into),colorarray(j))
        hold on
        plot(TOF_retentiont01l1,ones(1,length(TOF_retentiont01l1))*th(j),colorarray(j))
    end
    Height=max(th)*5;
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(:,1)),2*Height*ones(length(IntervalList01{ID_into}.intervallist_after7_combine(:,1)),1),'rd')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(:,2)),2*Height*ones(length(IntervalList01{ID_into}.intervallist_after7_combine(:,2)),1),'kd')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after6(:,1)),Height*ones(length(IntervalList01{ID_into}.intervallist_after6(:,1)),1),'ro')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after6(:,2)),Height*ones(length(IntervalList01{ID_into}.intervallist_after6(:,2)),1),'ko')
    grid on
    
    subplot(2,1,2)
    for j=1:6
        th(j)=TOF_getThreshold_range(TOF_XICs02{j}(:,ID_into),instrument_h_int, de_range,Times_noise_std);
        plot(TOF_retentiont02l1,TOF_XICs02{j}(:,ID_into),colorarray(j))
        hold on
        plot(TOF_retentiont02l1,ones(1,length(TOF_retentiont02l1))*th(j),colorarray(j))
    end
    Height=max(th)*5;
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(:,1)),2*Height*ones(length(IntervalList02{ID_into}.intervallist_after7_combine(:,1)),1),'rd')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(:,2)),2*Height*ones(length(IntervalList02{ID_into}.intervallist_after7_combine(:,2)),1),'kd')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after6(:,1)),Height*ones(length(IntervalList02{ID_into}.intervallist_after6(:,1)),1),'ro')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after6(:,2)),Height*ones(length(IntervalList02{ID_into}.intervallist_after6(:,2)),1),'ko')
    grid on
    
end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Generate all the AT AR AKL scores(not probability) for
%%%%%%%%%%%%%%%%%% only TOF data

for i=1:length(Pep_com)    
    iso_spec_pep=iso_com(i,:);    
    %%%%%%%%%%%%%%%%%%%% TOF to TOF contains AT AR AKL scores
    %%%%%%%%%%%%%%%%%%%%
    Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
    Interval_matrix02=IntervalList_com02{i}.intervallist_after7_combine;
    for k=1:size(Interval_matrix01,1)
        for j=1:6
            XICs01{k}(:,j)=TOF_XICs_com01{j}(Interval_matrix01(k,1):Interval_matrix01(k,2),i);
        end 
    end
    for k=1:size(Interval_matrix02,1)
        for j=1:6
            XICs02{k}(:,j)=TOF_XICs_com02{j}(Interval_matrix02(k,1):Interval_matrix02(k,2),i);
        end 
    end
%     for j=1:6
%         XICs01(:,j)=TOF_XICs_com01{j}(:,i);
%         XICs02(:,j)=TOF_XICs_com02{j}(:,i);
%     end
    %%%%%% intervals01 and intervals02 are matrix;
    %%%%%% XICs01 and XICs02 are matrixs (n*6) 6 colums are for different
    %%%%%% iso mz
    [TOFscorev1{i}.LOGAKL,TOFscorev1{i}.AR,TOFscorev1{i}.AT,TOFscorev1{i}.Normal_AT,TOFscorev1{i}.T01,TOFscorev1{i}.Normal_T01,TOFscorev1{i}.T02,TOFscorev1{i}.Normal_T02]=GenerateRTKL(Interval_matrix01,Interval_matrix02,XICs01,XICs02,TOF_retentiont01l1,TOF_retentiont02l1,iso_spec_pep);
    clear XICs01 XICs02
    %%%%%%%%%%%%%%%%%%%% 
end

save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\TOFscorev1 TOFscorev1
load D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\TOFscorev1
TOFscore=TOFscorev1;
save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\TOFscore TOFscore

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% generate AR AT AKL model based on the training data
%%%%%%%%%%%%%%%%%% 

%%%%%% 1. Choose AT based on AKL and AR
%%%%%% 2. Choose AKL based on AR and AT
%%%%%% 3. Choose AR based on AKL and AT

%%%%%%%%%%%%%%%%% Step 1 to get AT model 
Training_Row_Column_Num=[];
Timepair_corr=[];
Timepair_corr_normal=[];
Timepair_noncorr=[];
Timepair_noncorr_normal=[];
AT_corr=[];AT_non_corr=[];
Training_id_AT=[];
for i=1:length(TOFscore)
    
        V_LOGAKL=TOFscore{i}.LOGAKL;
        V_AR=TOFscore{i}.AR;
        V_Normal_AT=TOFscore{i}.Normal_AT;
 
        %%%%%%%%%% choose AT
        [AKL_R_good_id,AKL_C_good_id,AKL_V_good_id]=find(V_LOGAKL<=-6);
        [AR_R_good_id,AR_C_good_id,AR_V_good_id]=find(V_AR>=0.80);
             
        AKL_good_id_vector=[AKL_R_good_id,AKL_C_good_id];
        AR_good_id_vector=[AR_R_good_id,AR_C_good_id];
        Good_R_id=[];
        Good_C_id=[];
        if ~isempty(AKL_good_id_vector) && ~isempty(AR_good_id_vector)            
           
            for j=1:size(AKL_good_id_vector,1)
                index_diff01=AR_good_id_vector(:,1)-AKL_good_id_vector(j,1);
                index_diff02=AR_good_id_vector(:,2)-AKL_good_id_vector(j,2);
                Same_index=find(abs(index_diff01)+abs(index_diff02)==0);
                if ~isempty(Same_index)
                    Good_R_id=[Good_R_id;AKL_good_id_vector(j,1)];
                    Good_C_id=[Good_C_id;AKL_good_id_vector(j,2)];
                end
            end
            
            if ~isempty(Good_R_id)                
               
                R_corr=Good_R_id;
                C_corr=Good_C_id;
                
                if length(R_corr)==1 && length(C_corr)==1

                   if abs(V_Normal_AT(R_corr,C_corr))<=0.05                   
                    
                        Training_id_AT=[Training_id_AT; i];    
                        
                        Vector_T01_noncorr_Row=TOFscore{i}.T01(:,1);
                        Vector_T02_noncorr_Col=TOFscore{i}.T02(1,:);
                        Vector_T01_noncorr_Row(R_corr)=[];
                        Vector_T02_noncorr_Col(C_corr)=[];
                        Vector_T01_normal_noncorr_Row=TOFscore{i}.Normal_T01(:,1);
                        Vector_T02_normal_noncorr_Col=TOFscore{i}.Normal_T02(1,:);
                        Vector_T01_normal_noncorr_Row(R_corr)=[];
                        Vector_T02_normal_noncorr_Col(C_corr)=[];
                        Timepair_corr=[Timepair_corr;TOFscore{i}.T01(R_corr,C_corr),TOFscore{i}.T02(R_corr,C_corr)];
                        Timepair_corr_normal=[Timepair_corr_normal;TOFscore{i}.Normal_T01(R_corr,C_corr),TOFscore{i}.Normal_T02(R_corr,C_corr)];
                        Timepair_noncorr=[Timepair_noncorr;TOFscore{i}.T01(R_corr,C_corr)*ones(length(Vector_T02_noncorr_Col),1),Vector_T02_noncorr_Col'];
                        Timepair_noncorr_normal=[Timepair_noncorr_normal;TOFscore{i}.Normal_T01(R_corr,C_corr)*ones(length(Vector_T02_normal_noncorr_Col),1),Vector_T02_normal_noncorr_Col'];
                        
                        V_AT_corr=TOFscore{i}.Normal_AT(R_corr,C_corr);
                        Vector_AT_noncorr_Row=TOFscore{i}.Normal_AT(R_corr,:);
                        Vector_AT_noncorr_Col=TOFscore{i}.Normal_AT(:,C_corr);
                        Vector_AT_noncorr_Row(C_corr)=[];
                        Vector_AT_noncorr_Col(R_corr)=[];
                        AT_corr=[AT_corr;V_AT_corr];
                        AT_non_corr=[AT_non_corr;Vector_AT_noncorr_Row';Vector_AT_noncorr_Col];
                        Training_Row_Column_Num=[Training_Row_Column_Num; R_corr C_corr];
                   end                    
                    
                end
            end
        end
end
PP=polyfit(Timepair_corr_normal(:,1),Timepair_corr_normal(:,2),4);
figure
plot(Timepair_corr_normal(:,1),Timepair_corr_normal(:,2),'r.')
XXX=0:0.01:1;
YYY=polyval(PP,XXX);
hold on
plot(XXX,YYY)
grid on
AT_corr=polyval(PP,Timepair_corr_normal(:,1))-Timepair_corr_normal(:,2);
AT_non_corr=polyval(PP,Timepair_noncorr_normal(:,1))-Timepair_noncorr_normal(:,2);
[timemodel0102_MU_sig,timemodel0102_SIGMA_sig] = normfit(AT_corr);
[timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig] = normfit(AT_non_corr);
t=-1:0.001:1;
P_timenorm_sigv2=normpdf(t,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
P_timenorm_notsigv2=normpdf(t,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
xx=-1:0.001:1;
figure
subplot(2,2,1)
hist(AT_corr,xx);grid on;
subplot(2,2,3)
hist(AT_non_corr,xx);grid on;
subplot(1,2,2)
plot(t,P_timenorm_sigv2,'r');hold on
plot(t,P_timenorm_notsigv2)
grid on
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Step 2 to get AKL model
AKL0102_corr=[];AKL0102_non_corr=[];
Training_id_AKL=[];
for i=1:length(TOFscore)
    
        V_LOGAKL=TOFscore{i}.LOGAKL;
        V_AR=TOFscore{i}.AR;
        Score_V_Normal_AT_corr=normpdf(TOFscore{i}.Normal_AT,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
        Score_V_Normal_AT_non_corr=normpdf(TOFscore{i}.Normal_AT,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
        Score_AT=log(Score_V_Normal_AT_corr./Score_V_Normal_AT_non_corr);
        %%%%%%%%%% choose AKL
        [AT_R_good_id,AT_C_good_id,AT_V_good_id]=find(Score_AT>=0);
        [AR_R_good_id,AR_C_good_id,AR_V_good_id]=find(V_AR>=0.90);
             
        AT_good_id_vector=[AT_R_good_id,AT_C_good_id];
        AR_good_id_vector=[AR_R_good_id,AR_C_good_id];

        Good_R_id=[];
        Good_C_id=[];
        if ~isempty(AT_good_id_vector) && ~isempty(AR_good_id_vector)
            
           
            for j=1:size(AT_good_id_vector,1)
                index_diff01=AR_good_id_vector(:,1)-AT_good_id_vector(j,1);
                index_diff02=AR_good_id_vector(:,2)-AT_good_id_vector(j,2);
                Same_index=find(abs(index_diff01)+abs(index_diff02)==0);
                if ~isempty(Same_index)
                    Good_R_id=[Good_R_id;AT_good_id_vector(j,1)];
                    Good_C_id=[Good_C_id;AT_good_id_vector(j,2)];
                end
            end
            
            if ~isempty(Good_R_id)
                
               
                R_corr=Good_R_id;
                C_corr=Good_C_id;
                
                if length(R_corr)==1 && length(C_corr)==1
             
                        Training_id_AKL=[Training_id_AKL; i];    
                        
                        V_AKL_corr=TOFscore{i}.LOGAKL(R_corr,C_corr);
                        Vector_AKL_noncorr_Row=TOFscore{i}.LOGAKL(R_corr,:);
                        Vector_AKL_noncorr_Col=TOFscore{i}.LOGAKL(:,C_corr);
                        Vector_AKL_noncorr_Row(C_corr)=[];
                        Vector_AKL_noncorr_Col(R_corr)=[];
                        AKL0102_corr=[AKL0102_corr;V_AKL_corr];
                        AKL0102_non_corr=[AKL0102_non_corr;Vector_AKL_noncorr_Row';Vector_AKL_noncorr_Col];        
                
               end
            end
        end

end
[AKL0102_MU_sig,AKL0102_SIGMA_sig]=normfit(AKL0102_corr);
[AKL0102_MU_notsig,AKL0102_SIGMA_notsig]=normfit(AKL0102_non_corr);
xx=-20:0.1:10;
figure
subplot(2,2,1)
hist(AKL0102_corr,xx);grid on;
subplot(2,2,3)
hist(AKL0102_non_corr,xx);grid on;
subplot(1,2,2)
plot(xx,normpdf(xx,AKL0102_MU_sig,AKL0102_SIGMA_sig),'r')
hold on
plot(xx,normpdf(xx,AKL0102_MU_notsig,AKL0102_SIGMA_notsig));grid on;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Step 3 to get AR model
AR_corr=[];AR_non_corr=[];
Training_id_AR=[];
for i=1:length(TOFscore)
    
        V_LOGAKL=TOFscore{i}.LOGAKL;
        V_AR=TOFscore{i}.AR;
        Score_V_Normal_AT_corr=normpdf(TOFscore{i}.Normal_AT,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
        Score_V_Normal_AT_non_corr=normpdf(TOFscore{i}.Normal_AT,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
        Score_AT=log(Score_V_Normal_AT_corr./Score_V_Normal_AT_non_corr);
        %%%%%%%%%% choose AR
        [AT_R_good_id,AT_C_good_id,AT_V_good_id]=find(Score_AT>=0);
        [AKL_R_good_id,AKL_C_good_id,AKL_V_good_id]=find(V_LOGAKL<=-9);
             
        AT_good_id_vector=[AT_R_good_id,AT_C_good_id];
        AKL_good_id_vector=[AKL_R_good_id,AKL_C_good_id];

        Good_R_id=[];
        Good_C_id=[];
        if ~isempty(AT_good_id_vector) && ~isempty(AKL_good_id_vector)
            
           
            for j=1:size(AT_good_id_vector,1)
                index_diff01=AKL_good_id_vector(:,1)-AT_good_id_vector(j,1);
                index_diff02=AKL_good_id_vector(:,2)-AT_good_id_vector(j,2);
                Same_index=find(abs(index_diff01)+abs(index_diff02)==0);
                if ~isempty(Same_index)
                    Good_R_id=[Good_R_id;AT_good_id_vector(j,1)];
                    Good_C_id=[Good_C_id;AT_good_id_vector(j,2)];
                end
            end
            
            if ~isempty(Good_R_id)
                
               
                R_corr=Good_R_id;
                C_corr=Good_C_id;
                
                if length(R_corr)==1 && length(C_corr)==1
             
                        Training_id_AR=[Training_id_AR; i];    
                        
                    V_AR_corr=TOFscore{i}.AR(R_corr,C_corr);
                    Vector_AR_noncorr_Row=TOFscore{i}.AR(R_corr,:);
                    Vector_AR_noncorr_Col=TOFscore{i}.AR(:,C_corr);
                    Vector_AR_noncorr_Row(C_corr)=[];
                    Vector_AR_noncorr_Col(R_corr)=[];
                    AR_corr=[AR_corr;V_AR_corr];
                    AR_non_corr=[AR_non_corr;Vector_AR_noncorr_Row';Vector_AR_noncorr_Col];   
                    
               end
            end
        end

end
aa= AR_corr>=0 & AR_corr<=1;
AR_corr_noNaN=AR_corr(aa);
aa= AR_non_corr>=0 & AR_non_corr<=1;
AR_non_corr_noNaN=AR_non_corr(aa);
aa= AR_non_corr>=0 & AR_non_corr<=1;
AR_non_corr_noNaN=AR_non_corr(aa);
AR0102_PARMHAT1=betafit(AR_corr_noNaN);
AR0102_PARMHAT2=betafit(AR_non_corr_noNaN);
b1=0:0.05:1;
Y1_sig= betapdf(b1,AR0102_PARMHAT1(1),AR0102_PARMHAT1(2));
b1=0:0.05:1;
Y1_notsig= betapdf(b1,AR0102_PARMHAT2(1),AR0102_PARMHAT2(2));
figure
plot(b1,Y1_sig,'r');hold on;plot(b1,Y1_notsig)
%%%%%%%%%%%%%%%%%%

length(Training_id_AKL)
length(Training_id_AR)
length(Training_id_AT)

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% check training data
%%%%%%%%%%%%%%%
clear Interval_data01
clear log_KL_Value_corr
clear Interval_data01 Interval_data02
clear log_KL_Value_noncorr01
num=1;
Training_id=Training_id_AT;
for i=1:length(Training_id)
    data01_intervalid=Training_Row_Column_Num(i,1);
    data02_intervalid=Training_Row_Column_Num(i,2);
    ID_into=Com_ID_0102(Training_id(i));
    scan_start01=IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,1);
    scan_end01=IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,2);
    for j=1:2
        Interval_data01(:,j)=TOF_XICs01{j}(scan_start01:scan_end01,ID_into);
    end
    scan_start02=IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,1);
    scan_end02=IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,2);
    for j=1:2
        Interval_data02(:,j)=TOF_XICs02{j}(scan_start02:scan_end02,ID_into);
    end
    
    Interval_Total=[Interval_data01;Interval_data02];
%     Interval_Total=Interval_data01;
    Weight_ISO=sum(Interval_Total,2)/sum(sum(Interval_Total));
    TOTAL=sum(diag(Weight_ISO)*Interval_Total,1);
    Obe_TOFISO_Pattern=TOTAL/sum(TOTAL);
    
    Th_ISO_pattern=iso(ID_into,1:2)/sum(iso(ID_into,1:2));
    log_KL_Value_corr(i)=KL_calculate(Obe_TOFISO_Pattern,Th_ISO_pattern);
    
    clear Interval_data01 Interval_data02
    
%     ID01=1:size(IntervalList01{ID_into}.intervallist_after7_combine,1);
    if size(IntervalList01{ID_into}.intervallist_after7_combine,1)~=1
        if data01_intervalid==1;
            ID01=2;
        else if data01_intervalid==size(IntervalList01{ID_into}.intervallist_after7_combine,1)
                ID01=size(IntervalList01{ID_into}.intervallist_after7_combine,1)-1;
            else
                ID01=[data01_intervalid-1;data01_intervalid+1];
            end
        end
%         ID01(data01_intervalid)=[];
        for k=1:length(ID01)
            scan_start01=IntervalList01{ID_into}.intervallist_after7_combine(ID01(k),1);
            scan_end01=IntervalList01{ID_into}.intervallist_after7_combine(ID01(k),2);
            for j=1:2
                Interval_data01(:,j)=TOF_XICs01{j}(scan_start01:scan_end01,ID_into);
            end
            Weight_ISO=sum(Interval_data01,2)/sum(sum(Interval_data01));
            TOTAL=sum(diag(Weight_ISO)*Interval_data01,1);
            Obe_TOFISO_Pattern=TOTAL/sum(TOTAL);

            Th_ISO_pattern=iso(ID_into,1:2)/sum(iso(ID_into,1:2));
            log_KL_Value_noncorr01(num)=KL_calculate(Obe_TOFISO_Pattern,Th_ISO_pattern);
            num=num+1;
            clear Interval_data01
        end
        clear Interval_data01
    end
end
figure;hist(log_KL_Value_corr)
figure;hist(log_KL_Value_noncorr01)
var(log_KL_Value_corr)

Times_noise_std=3;%6;
instrument_h_int=10^7;
de_range=10^5;%10^4;
colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
for i=36:40%%length(posi01v2)
    data01_intervalid=Training_Row_Column_Num(i,1);
    data02_intervalid=Training_Row_Column_Num(i,2);
    ID_into=Com_ID_0102(Training_id(i));
    
    figure
    subplot(2,1,1)
    for j=1:6
        th(j)=TOF_getThreshold_range(TOF_XICs01{j}(:,ID_into),instrument_h_int, de_range,Times_noise_std);
        plot(TOF_retentiont01l1,TOF_XICs01{j}(:,ID_into),colorarray(j))
        hold on
        plot(TOF_retentiont01l1,ones(1,length(TOF_retentiont01l1))*th(j),colorarray(j))
    end
    Height=max(th)*5;
    
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,1)),3*Height*ones(length(IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,1)),1),'r*')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,2)),3*Height*ones(length(IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,2)),1),'k*')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(:,1)),2*Height*ones(length(IntervalList01{ID_into}.intervallist_after7_combine(:,1)),1),'rd')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(:,2)),2*Height*ones(length(IntervalList01{ID_into}.intervallist_after7_combine(:,2)),1),'kd')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after6(:,1)),Height*ones(length(IntervalList01{ID_into}.intervallist_after6(:,1)),1),'ro')
    stem(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after6(:,2)),Height*ones(length(IntervalList01{ID_into}.intervallist_after6(:,2)),1),'ko')
    grid on
    
    subplot(2,1,2)
    for j=1:6
        th(j)=TOF_getThreshold_range(TOF_XICs02{j}(:,ID_into),instrument_h_int, de_range,Times_noise_std);
        plot(TOF_retentiont02l1,TOF_XICs02{j}(:,ID_into),colorarray(j))
        hold on
        plot(TOF_retentiont02l1,ones(1,length(TOF_retentiont02l1))*th(j),colorarray(j))
    end
    Height=max(th)*5;
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,1)),3*Height*ones(length(IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,1)),1),'r*')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,2)),3*Height*ones(length(IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,2)),1),'k*')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(:,1)),2*Height*ones(length(IntervalList02{ID_into}.intervallist_after7_combine(:,1)),1),'rd')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(:,2)),2*Height*ones(length(IntervalList02{ID_into}.intervallist_after7_combine(:,2)),1),'kd')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after6(:,1)),Height*ones(length(IntervalList02{ID_into}.intervallist_after6(:,1)),1),'ro')
    stem(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after6(:,2)),Height*ones(length(IntervalList02{ID_into}.intervallist_after6(:,2)),1),'ko')
    grid on    

    
end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Generate Testing data
%%%%%%%%%%%%%%%%%% The Testing data is the whole data for aligning

Testing_id=1:length(TOFscorev1);
% Testing_id(Training_id)=[];
Testingdata=TOFscorev1(Testing_id);

% iso_training=iso_com(Training_id);
% Pep_training=Pep_com(Training_id);
% TOF_mass_training=TOF_mass_com(Training_id);
% TOF_charge_training=TOF_charge_com(Training_id);
% Protein_training=Protein_matrix_final_com(Training_id);
% for i=1:6
%     TOF_XICs_training01{i}=TOF_XICs_com01{i}(:,Training_id);
% end
% IntervalList_training01=IntervalList_com01(Training_id);
% for i=1:6
%     TOF_XICs_training02{i}=TOF_XICs_com02{i}(:,Training_id);
% end
% IntervalList_training02=IntervalList_com02(Training_id);

iso_testing=iso_com(Testing_id,:);
Pep_testing=Pep_com(Testing_id);
TOF_mass_testing=TOF_mass_com(Testing_id);
TOF_charge_testing=TOF_charge_com(Testing_id);
Protein_testing=Protein_matrix_final_com(Testing_id);
for i=1:6
    TOF_XICs_testing01{i}=TOF_XICs_com01{i}(:,Testing_id);
end
IntervalList_testing01=IntervalList_com01(Testing_id);
for i=1:6
    TOF_XICs_testing02{i}=TOF_XICs_com02{i}(:,Testing_id);
end
IntervalList_testing02=IntervalList_com02(Testing_id);
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% corr>non_corr
num=0;
num_notdetect=0;
TOF_goodalign_corrBnon_on_testing=[];

for i=1:length(Testingdata) 
        AR_corr_Prob=betapdf(Testingdata{i}.AR,AR0102_PARMHAT1(1),AR0102_PARMHAT1(2));
        AR_non_corr_Prob=betapdf(Testingdata{i}.AR,AR0102_PARMHAT2(1),AR0102_PARMHAT2(2));
        AKL_corr_Prob=normpdf(Testingdata{i}.LOGAKL,AKL0102_MU_sig,AKL0102_SIGMA_sig);
        AKL_non_corr_Prob=normpdf(Testingdata{i}.LOGAKL,AKL0102_MU_notsig,AKL0102_SIGMA_notsig);
        AT_corr_Prob=normpdf(polyval(PP,Testingdata{i}.Normal_T01)-Testingdata{i}.Normal_T02,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
        AT_non_corr_Prob=normpdf(polyval(PP,Testingdata{i}.Normal_T01)-Testingdata{i}.Normal_T02,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
        Prob_AR_matrix=log(AR_corr_Prob./AR_non_corr_Prob);
        Prob_AKL_matrix=log(AKL_corr_Prob./AKL_non_corr_Prob);
        Prob_AT_matrix=log(AT_corr_Prob./AT_non_corr_Prob);
%         Total_Prob_matrix=Prob_AR_matrix+Prob_AKL_matrix+Prob_AT_matrix;
%         Total_Prob_matrix=Prob_AKL_matrix+Prob_AT_matrix;
        Total_Prob_matrix=Prob_AT_matrix;
%         Total_Prob_matrix=Prob_AKL_matrix;
        R_detect_row=[];
        C_detect_row=[];       
        for R=1:size(Total_Prob_matrix,1)
            [V_m,I_m]=max(Total_Prob_matrix(R,:));
            if V_m>=0
                R_detect_row=[R_detect_row;R];
                C_detect_row=[C_detect_row;I_m];
            end
        end
        show=[R_detect_row,C_detect_row];
        R_detect=[];
        C_detect=[];
        for C=1:size(Total_Prob_matrix,2)            
            ID_C=find(C_detect_row==C);
            if ~isempty(ID_C)
                [V_m,I_m]=max(Total_Prob_matrix(R_detect_row(ID_C),C));
                if V_m>=0
                    R_detect=[R_detect;R_detect_row(ID_C(I_m))];
                    C_detect=[C_detect;C_detect_row(ID_C(I_m))];                    
                end
            end        
        end
        show=[R_detect,C_detect];
%         [R_detect,C_detect,V_detect]=find(Total_Prob_matrix>=0);
        [V_max_detect,I_max_detect]=max(Total_Prob_matrix);
        [V_max_detectv1,I_max_detectv1]=max(V_max_detect);    
        R_max_detect=I_max_detect(I_max_detectv1);
        C_max_detect=I_max_detectv1;
        if ~isempty(R_detect)
            num=num+1;
            TOF_goodalign_corrBnon_on_testing=[TOF_goodalign_corrBnon_on_testing;i];
            Pepinformdetect(i).ID_detect=i;
            Pepinformdetect(i).R_detect=R_detect;
            Pepinformdetect(i).C_detect=C_detect;
            Pepinformdetect(i).R_max_detect=R_max_detect;
            Pepinformdetect(i).C_max_detect=C_max_detect;
        else
            Pepinformdetect(i).ID_detect=0;
            Pepinformdetect(i).R_detect=0;
            Pepinformdetect(i).C_detect=0;
            Pepinformdetect(i).R_max_detect=0;
            Pepinformdetect(i).C_max_detect=0;            
        end
     
end
length(unique(TOF_goodalign_corrBnon_on_testing))

% save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\Testingdata Testingdata
% save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\Pepinformdetect Pepinformdetect
% save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\Pepinformdetect_AT Pepinformdetect
% save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\Pepinformdetect_AKL Pepinformdetect
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% check the alignment result
Times_noise_std=3;%6;
instrument_h_int=10^7;
de_range=10^5;%10^4;
colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
for i=101:105%length(Pepinformdetect)
    
    Spec_Testing_id=Pepinformdetect(i).ID_detect;
    Detect_interval_Id01=Pepinformdetect(i).R_detect;
    Detect_interval_Id02=Pepinformdetect(i).C_detect;
    clear Interval_data01 Interval_data02
    for j=1:6
        Interval_data01(:,j)=TOF_XICs_testing01{j}(:,Spec_Testing_id);
        Interval_data02(:,j)=TOF_XICs_testing02{j}(:,Spec_Testing_id);        
    end
    figure
    subplot(2,1,1)
    for j=1:6
        th(j)=TOF_getThreshold_range(Interval_data01(:,j),instrument_h_int, de_range,Times_noise_std);
        plot(TOF_retentiont01l1,Interval_data01(:,j),colorarray(j))
        hold on
        plot(TOF_retentiont01l1,ones(1,length(TOF_retentiont01l1))*th(j),colorarray(j))
    end
    Height=max(th)*5;    
    stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id01,1)),4*Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id01,1)),1),'r*')
    stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id01,2)),4*Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id01,2)),1),'k*')
%     stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(:,1)),3*Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(:,1)),1),'rd')
%     stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(:,2)),3*Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(:,2)),1),'kd')
%     stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7(:,1)),2*Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after7(:,1)),1),'r+')
%     stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7(:,2)),2*Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after7(:,2)),1),'k+')
%     stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after6(:,1)),Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after6(:,1)),1),'ro')
%     stem(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after6(:,2)),Height*ones(length(IntervalList_testing01{Spec_Testing_id}.intervallist_after6(:,2)),1),'ko')

    for kkk=1:length(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id01,1))
        interval_mid_time=(TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id01(kkk),1))+TOF_retentiont01l1(IntervalList_testing01{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id01(kkk),2)))/2;
        text(interval_mid_time,3*Height,num2str(kkk));
    end    
    grid on
    
    subplot(2,1,2)
    for j=1:6
        th(j)=TOF_getThreshold_range(Interval_data02(:,j),instrument_h_int, de_range,Times_noise_std);
        plot(TOF_retentiont02l1,Interval_data02(:,j),colorarray(j))
        hold on
        plot(TOF_retentiont02l1,ones(1,length(TOF_retentiont02l1))*th(j),colorarray(j))
    end
    Height=max(th)*5;    
    stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id02,1)),4*Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id02,1)),1),'r*')
    stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id02,2)),4*Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id02,2)),1),'k*')
%     stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(:,1)),3*Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(:,1)),1),'rd')
%     stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(:,2)),3*Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(:,2)),1),'kd')
%     stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7(:,1)),2*Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after7(:,1)),1),'r+')
%     stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7(:,2)),2*Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after7(:,2)),1),'k+')
%     stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after6(:,1)),Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after6(:,1)),1),'ro')
%     stem(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after6(:,2)),Height*ones(length(IntervalList_testing02{Spec_Testing_id}.intervallist_after6(:,2)),1),'ko')
%   
    
    for kkk=1:length(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id02,1))
        interval_mid_time=(TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id02(kkk),1))+TOF_retentiont02l1(IntervalList_testing02{Spec_Testing_id}.intervallist_after7_combine(Detect_interval_Id02(kkk),2)))/2;
        text(interval_mid_time,3*Height,num2str(kkk));
    end    
    grid on
    clear Interval_data01 Interval_data02
end
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Generate Align Result structure
for i=1:length(Pepinformdetect)
    TOF_align_Result{i}.peptideseqence=Pep_testing{i};
    TOF_align_Result{i}.mass=TOF_mass_testing(i);
    TOF_align_Result{i}.chargestate=TOF_charge_testing(i);
    TOF_align_Result{i}.iso=iso_testing(i,:);
    TOF_align_Result{i}.Protein=Protein_testing{i};
    
    if Pepinformdetect(i).ID_detect==0
        TOF_align_Result{i}.Judge_Align=0;        
    else
        TOF_align_Result{i}.Judge_Align=1;  
        Interval_id_A=Pepinformdetect(i).R_detect;
        for j=1:length(Interval_id_A)
            Interval=IntervalList_testing01{i}.intervallist_after7_combine(Interval_id_A(j),:);
            TOF_align_Result{i}.Intervallist_A(j,:)=Interval;
            TOF_align_Result{i}.TimeInterval_A(j,:)=TOF_retentiont01l1(Interval);
            TOF_align_Result{i}.alignindex_A(j,1)=j;
            for k=1:6
                TOF_align_Result{i}.Interval_data_A{j}(:,k)=TOF_XICs_testing01{k}(Interval(1):Interval(2),i);
            end
        end

        Interval_id_B=Pepinformdetect(i).C_detect;
        for j=1:length(Interval_id_B)
            Interval=IntervalList_testing02{i}.intervallist_after7_combine(Interval_id_B(j),:);
            TOF_align_Result{i}.Intervallist_B(j,:)=Interval;
            TOF_align_Result{i}.TimeInterval_B(j,:)=TOF_retentiont02l1(Interval);
            TOF_align_Result{i}.alignindex_B(j,1)=j;
            for k=1:6
                TOF_align_Result{i}.Interval_data_B{j}(:,k)=TOF_XICs_testing02{k}(Interval(1):Interval(2),i);
            end
        end
    end
end
% save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\TOF_align_Result TOF_align_Result
% save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\TOF_align_Result_AT TOF_align_Result
% save D:\Program\QTOF_replicate_identification\TOF_TOF_MATfile\TOF_align_Result_AKL TOF_align_Result
%%%%%%%%%%%%%%%%%%%
num=1;
for i=1:length(TOF_align_Result)
    if TOF_align_Result{i}.Judge_Align==1
        for k=1:length(TOF_align_Result{i}.Interval_data_A)
            TestingResult01.intervalsdata=TOF_align_Result{i}.Interval_data_A{k};
            TestingResult01.iso=TOF_align_Result{i}.iso;
            TestingResult02.intervalsdata=TOF_align_Result{i}.Interval_data_B{k};
            TestingResult02.iso=TOF_align_Result{i}.iso;
            [finalO18rate01(num,:),finalO18f01(num,:)]=OrbitrapProduceO18RatesV1(TestingResult01); 
            [finalO18rate02(num,:),finalO18f02(num,:)]=OrbitrapProduceO18RatesV1(TestingResult02);
            num=num+1;
        end   
    end
end
std(log2(finalO18rate01(:,3)))
std(log2(finalO18rate02(:,3)))
std(log2(finalO18rate01(:,3))-log2(finalO18rate02(:,3)))
%%%%%%%%%%%%%%%%%%%






