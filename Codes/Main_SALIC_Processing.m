%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This program is for processing SILAC SC10 data. Finish the interval
%%%% detection and calculate the L/H ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

Maxquant_file_path='C:\Users\Jian Cui\Desktop\Haskin_raw_data10\peptides.txt';%% SC25
% Maxquant_file_path='C:\Users\Jian Cui\Desktop\Xiaolin Expt01212012\A+B10_Maxquant\combined\txt\peptides.txt';%% AB10
[Maxquant_SC10data, Maxquant_SC10result]=readtext(Maxquant_file_path,'\t');
Maxquant_Razor_Protein=Maxquant_SC10data(2:end,34);
Maxquant_peptide_sequence=Maxquant_SC10data(2:end,8);
Maxquant_Razor_Protein_HLR=Maxquant_SC10data(2:end,43);
Maxquant_HLRatio=[];
for i=1:length(Maxquant_Razor_Protein_HLR)
    if ~isnan(Maxquant_Razor_Protein_HLR{i})
        Maxquant_HLRatio=[Maxquant_HLRatio; Maxquant_Razor_Protein_HLR{i},i];        
    end    
end
size(Maxquant_HLRatio)
mean(Maxquant_HLRatio(:,1))
length(unique(Maxquant_peptide_sequence))
figure
hist(log2(Maxquant_HLRatio(:,1)),[-20:0.1:10]);grid on

num=1;
for i=1:length(Maxquant_Razor_Protein_HLR)
    if ~isnan(Maxquant_Razor_Protein_HLR{i})        
        Final_Maxquant_Razor_Protein{num}=Maxquant_Razor_Protein{i};
        Final_Maxquant_Razor_Protein_HLR(num)=Maxquant_Razor_Protein_HLR{i};        
        num=num+1;
    end
end
Final_Maxquant_Razor_Protein_unique=unique(Final_Maxquant_Razor_Protein);
for i=1:length(Final_Maxquant_Razor_Protein_unique)
    ID=strcmp(Final_Maxquant_Razor_Protein_unique{i},Final_Maxquant_Razor_Protein);
    Final_Maxquant_Razor_Protein_unique_HLR(i)=sum(Final_Maxquant_Razor_Protein_HLR(ID))/sum(ID);
end
length(Final_Maxquant_Razor_Protein_unique_HLR)
mean(log2(Final_Maxquant_Razor_Protein_unique_HLR))
std(log2(Final_Maxquant_Razor_Protein_unique_HLR))
figure
hist(log2(Final_Maxquant_Razor_Protein_unique_HLR),[-20:0.1:10])

clear Final_Maxquant_Razor_Protein Final_Maxquant_Razor_Protein_HLR

load D:\Program\QTOF_replicate_identification\haskin_new_salic_data\Protein_SC10_information
LHR_Prot=[];
for i=1:length(Protein_SC10_information.LHR)
    LHR_Prot=[LHR_Prot,Protein_SC10_information.LHR(i)];
    Protein_SC10_seq{i}=Protein_SC10_information.Accession{i};
    Protein_SC10_hitnum(i)=Protein_SC10_information.protHit(i);
end
ID_LHR_Prot=~isnan(LHR_Prot);
mean(log2(1./LHR_Prot(ID_LHR_Prot)))
std(log2(1./LHR_Prot(ID_LHR_Prot)))
figure
hist(log2(1./LHR_Prot(ID_LHR_Prot)),[-20:0.1:10])
length(LHR_Prot(ID_LHR_Prot))
Protein_SC10_LHR=LHR_Prot(ID_LHR_Prot);
Protein_SC10_seq_withLHR=Protein_SC10_seq(ID_LHR_Prot);
Protein_SC10_hitnum=Protein_SC10_hitnum(ID_LHR_Prot);

unique(Protein_SC10_seq_withLHR)
%%%%%%%%%%%%%

[SC10data,SC10result]=readtext('D:\Program\QTOF_replicate_identification\haskin_new_salic_data\sc10.txt','\t');

Peptide_SC10_information.prot_hit_num=SC10data(2:end,1);
Peptide_SC10_information.prot_acc=SC10data(2:end,2);
Peptide_SC10_information.prot_desc=SC10data(2:end,3);
Peptide_SC10_information.prot_score=SC10data(2:end,4);
Peptide_SC10_information.prot_mass=SC10data(2:end,5);
Peptide_SC10_information.prot_matches=SC10data(2:end,6);
Peptide_SC10_information.prot_matches_sig=SC10data(2:end,7);
Peptide_SC10_information.prot_sequences=SC10data(2:end,8);
Peptide_SC10_information.prot_sequences_sig=SC10data(2:end,9);
Peptide_SC10_information.prot_seq=SC10data(2:end,15);
Peptide_SC10_information.pep_isunique=SC10data(2:end,19);
Peptide_SC10_information.pep_exp_mz=SC10data(2:end,20);
Peptide_SC10_information.pep_exp_mr=SC10data(2:end,21);
Peptide_SC10_information.pep_exp_z=SC10data(2:end,22);
Peptide_SC10_information.pep_calc_mr=SC10data(2:end,23);
Peptide_SC10_information.pep_start=SC10data(2:end,25);
Peptide_SC10_information.pep_end=SC10data(2:end,26);
Peptide_SC10_information.pep_res_before=SC10data(2:end,32);
Peptide_SC10_information.pep_seq=SC10data(2:end,33);
Peptide_SC10_information.pep_res_after=SC10data(2:end,34);
Peptide_SC10_information.pep_num_match=SC10data(2:end,38);
Peptide_SC10_information.pep_scan_title=SC10data(2:end,39);
Peptide_SC10_information.pep_var_mod=SC10data(2:end,36);

for i=1:length(Peptide_SC10_information.pep_scan_title)
    
    Time_Text_Information=Peptide_SC10_information.pep_scan_title{i};
    ID_S=find(Time_Text_Information=='[');    
    Time_Text_Information(ID_S:end)=[];
    ID_ST=find(Time_Text_Information=='(');
    ID_ED=find(Time_Text_Information==')');
    for k=1:length(ID_ST)
        ID_equ=find(Time_Text_Information(ID_ST(k):ID_ED(k))=='=');
        Time(k)=str2num(Time_Text_Information(ID_ST(k)+ID_equ:ID_ED(k)-1));     
        ID_W=find(Time_Text_Information(1:ID_ST(k)-1)>='A');
        SEQ_cand=Time_Text_Information(ID_W(end):ID_ST(k)-1);
        Scan(k)=str2num(SEQ_cand(regexp(SEQ_cand,'\d')));
    end
    Peptide_SC10_information.pep_scan_information(i).Time=Time*60;
    Peptide_SC10_information.pep_scan_information(i).Scan=Scan;
    
    clear Time Scan    

end
save D:\Program\QTOF_replicate_identification\haskin_new_salic_data\Peptide_SC10_information Peptide_SC10_information

load D:\Program\QTOF_replicate_identification\haskin_new_salic_data\Peptide_SC10_information


[SC10retentiontl1,SC10datal1,SC10peakl,SC10retentiont,SC10MZInt_l1l2]=readrawdata('D:\Program\QTOF_replicate_identification\haskin_new_salic_data\SC10.mzXML');
save D:\Program\QTOF_replicate_identification\haskin_new_salic_data\SC10_information SC10retentiontl1 SC10datal1 SC10peakl SC10retentiont SC10MZInt_l1l2

load D:\Program\QTOF_replicate_identification\haskin_new_salic_data\SC10_information

for i=1:length(Peptide_SC10_information.pep_seq)
    SEQ_B=Peptide_SC10_information.pep_res_before{i};
    SEQ_M=Peptide_SC10_information.pep_seq{i};
    SEQ_A=Peptide_SC10_information.pep_res_after{i};
    SEQ=[SEQ_B,'.',SEQ_M,'.',SEQ_A];
    [peptidenew, massdiffList, isHeavy]=modprocess({SEQ});
    peptidenew=peptidenew{1};
    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,4,getmodificationformula(SEQ));
    iso(i,:)=isotopepattern;
%     csvalue=Peptide_SC10_information.pep_exp_z{i};
%     mzvalue=(weight+csvalue*1.0073)/csvalue;
%     Peptide_SC10_information.pep_calc_mr{i}
%     weight
end

SALIC_totalmzList=[];
Heavy_ID=[];
for i=1:length(Peptide_SC10_information.pep_seq)
    
    Pep_seq=Peptide_SC10_information.pep_seq{i};
    K_posi=find(Pep_seq=='K');
    R_posi=find(Pep_seq=='R');
    Mod_inform=Peptide_SC10_information.pep_var_mod{i};
    K_Mod_string='Label:13C(6) (K)';
    R_Mod_string='Label:13C(6) (R)';
    ID_K_Mod_str_match=strfind(Mod_inform,K_Mod_string);
    ID_R_Mod_str_match=strfind(Mod_inform,R_Mod_string);    
    N_K=length(K_posi);
    N_R=length(R_posi);
    if ~isempty(ID_K_Mod_str_match)        
        if ID_K_Mod_str_match-2>0
            Num_mod=Mod_inform(ID_K_Mod_str_match-2);
            if ~isempty(str2num(Num_mod))
                N_now_Kmod=str2num(Num_mod);
            else
                N_now_Kmod=1;
            end
        else
            N_now_Kmod=1;
        end
    else
        N_now_Kmod=0;
    end
    if ~isempty(ID_R_Mod_str_match)
        if ID_R_Mod_str_match-2>0
            Num_mod=Mod_inform(ID_R_Mod_str_match-2);
            if ~isempty(str2num(Num_mod))
                N_now_Rmod=str2num(Num_mod);
            else
                N_now_Rmod=1;
            end
        else
            N_now_Rmod=1;
        end
    else
        N_now_Rmod=0;
    end
    if ~isempty(ID_R_Mod_str_match) || ~isempty(ID_K_Mod_str_match)
        Heavy_ID=[Heavy_ID;i];
    end
    
%     N_K
%     N_R
%     N_now_Kmod
%     N_now_Rmod
        
    weight_now=Peptide_SC10_information.pep_calc_mr{i};
    csvalue=Peptide_SC10_information.pep_exp_z{i};

    weight=weight_now-(N_now_Kmod+N_now_Rmod)*6.020129;    
    mzvalue=(weight+csvalue*1.0073)/csvalue;
    SALIC_totalmzList=[SALIC_totalmzList; mzvalue mzvalue+(13.0034-12)/csvalue mzvalue+1.0034*2/csvalue mzvalue+1.0034*3/csvalue mzvalue+(N_K+N_R)*6.020129/csvalue mzvalue+(13.0034-12)/csvalue+(N_K+N_R)*6.020129/csvalue mzvalue+1.0034*2/csvalue+(N_K+N_R)*6.020129/csvalue  mzvalue+1.0034*3/csvalue+(N_K+N_R)*6.020129/csvalue];
end

tolerance=20;
for i=1:8
    SILAC_SC10_XICs_Tolerance{i}=getXIC_LC_new(SC10peakl,SALIC_totalmzList(:,i),tolerance);
end
save D:\Program\QTOF_replicate_identification\haskin_new_salic_data\SILAC_SC10_XICs_Tolerance SILAC_SC10_XICs_Tolerance
load D:\Program\QTOF_replicate_identification\haskin_new_salic_data\SILAC_SC10_XICs_Tolerance
colorarray=['r', 'k', 'g', 'b', 'm', 'y','c','g--'];
height=1000000;
for kk=251:255
    if length(Peptide_SC10_information.pep_scan_information(kk).Time)>1
        Scan_time=sum(Peptide_SC10_information.pep_scan_information(kk).Time)/2;
    else
        Scan_time=Peptide_SC10_information.pep_scan_information(kk).Time;
    end
%     Scan_end=Peptide_SC10_information.pep_scan_information(kk).Time;
    figure
    for i=1:8
        plot(SC10retentiontl1,SILAC_SC10_XICs_Tolerance{i}(:,kk),colorarray(i))
        hold on
    end
    stem(Scan_time,height,'ro');
%     stem(SC10retentiontl1(Scan_end),height,'ko');
    grid on
end
%%%%%%%%%%%%% SILAC interval detection

for i=1:length(Peptide_SC10_information.pep_seq)    
    pep01{i}=[Peptide_SC10_information.pep_res_before{i},'.',Peptide_SC10_information.pep_seq{i},'.',Peptide_SC10_information.pep_res_after{i}];
    Pep_unique(i)=Peptide_SC10_information.pep_isunique{i};
end
ID_unique=find(Pep_unique==1);
length(unique(pep01(ID_unique)))

posi01v1=[];
posi01v2=[];
posi01=[];
for i=1:length(pep01)
    SILAC_XICs01_sepc_pep=zeros(size(SILAC_SC10_XICs_Tolerance{1},1),8);
    for j=1:8
        SILAC_XICs01_sepc_pep(:,j)=SILAC_SC10_XICs_Tolerance{j}(:,i);        
    end
    [IntervalList,Jud_detect_good]=SILAC_Verification_intervaldetection_Spec_Pep(pep01{i},iso(i,:),SILAC_XICs01_sepc_pep);
    IntervalList01{i}=IntervalList;
    if Jud_detect_good==1
        posi01v1=[posi01v1;i];
%         if IntervalList.intervallist_after7_combine(1,1)~=0 && IntervalList.intervallist_after7_combine(1,2)~=0
%             posi01v2=[posi01v2;i];
%             SC10retentiontl1v1=SC10retentiontl1';
%             T_m=SC10retentiontl1v1(IntervalList.intervallist_after7_combine)-sum(Peptide_SC10_information.pep_scan_information(i).Time)/length(Peptide_SC10_information.pep_scan_information(i).Time);
%             T_v=T_m(:,1).*T_m(:,2);
%             Id_ms2interval=find(T_v<=0);
%             if ~isempty(Id_ms2interval)
%                 posi01=[posi01;i];
%             end            
%         end
    end    
end
posi01v1=[];
posi01v2=[];
posi01=[];
for i=1:length(IntervalList01)
    if IntervalList01{i}.intervallist_after7_combine(1,1)~=0 && IntervalList01{i}.intervallist_after7_combine(1,2)~=0
        posi01v2=[posi01v2;i];
        SC10retentiontl1v1=SC10retentiontl1';
        T_m=SC10retentiontl1v1(IntervalList01{i}.intervallist_after7_combine)-sum(Peptide_SC10_information.pep_scan_information(i).Time)/length(Peptide_SC10_information.pep_scan_information(i).Time);
        T_v=T_m(:,1).*T_m(:,2);
        Id_ms2interval=find(T_v<=0);
        if ~isempty(Id_ms2interval)
            posi01=[posi01;i Id_ms2interval];
        end            
    end
end

save D:\Program\QTOF_replicate_identification\haskin_new_salic_data\IntervalList01_vnew IntervalList01
load D:\Program\QTOF_replicate_identification\haskin_new_salic_data\IntervalList01_vnew 
%%%%%%%%%%%%%%%%%%%%%%%%% check
colorarray=['r', 'k', 'g', 'b', 'm', 'y','c','r:'];
height=1000000;
for kkk=204:206%1:length(posi01)
    kk=posi01(kkk,1);
    pep_seq_string=Peptide_SC10_information.pep_seq{kk};
    if length(Peptide_SC10_information.pep_scan_information(kk).Time)>1
        Scan_time=sum(Peptide_SC10_information.pep_scan_information(kk).Time)/2;
    else
        Scan_time=Peptide_SC10_information.pep_scan_information(kk).Time;
    end
%     Scan_end=Peptide_SC10_information.pep_scan_information(kk).Time;
    figure
    for i=1:8
        plot(SC10retentiontl1,SILAC_SC10_XICs_Tolerance{i}(:,kk),colorarray(i))
        hold on
    end
    stem(Scan_time,height,'ro');
    grid on
    
    stem(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1)),height*2*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'r*')
    stem(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)),height*2*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'k*')
    X=(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1))+SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)))/2;
    for RRR=1:length(X)
        text(X(RRR),height*3,num2str(IntervalList01{kk}.Labelling_efficiency_after7_combine(RRR)))
    end
    title(pep_seq_string)

end

HLRATIO=zeros(1,length(posi01));
Protein_unique_id=zeros(1,length(posi01));
for kkk=1:length(posi01)
    kk=posi01(kkk,1);
    ID_c=posi01(kkk,2);    
    Protein_unique_id(kkk)=Peptide_SC10_information.pep_isunique{kk};    
    Protein_hitnum(kkk)=Peptide_SC10_information.prot_hit_num{kk};
    Peptide_seq_withmod{kkk}=[Peptide_SC10_information.pep_seq{kk},Peptide_SC10_information.pep_var_mod{kk}];
    Peptide_seq{kkk}=Peptide_SC10_information.pep_seq{kk};
    Protein_seq{kkk}=Peptide_SC10_information.prot_acc{kk};
    HLRATIO(kkk)=IntervalList01{kk}.Labelling_efficiency_after7_combine(ID_c); 
    IntervalScan=IntervalList01{kk}.intervallist_after7_combine(ID_c,:);
    SILAC_XICs01_sepc_pep=zeros(size(SILAC_SC10_XICs_Tolerance{1},1),8);
    for j=1:8
        SILAC_XICs01_sepc_pep(:,j)=SILAC_SC10_XICs_Tolerance{j}(:,kk);        
    end
    figure
    for j=1:8
        plot(SILAC_XICs01_sepc_pep(IntervalScan(1):IntervalScan(2),j),colorarray(j))
        hold on
    end
    grid on
    
    LightMatrix=SILAC_XICs01_sepc_pep(IntervalScan(1):IntervalScan(2),1:4);
    HeavyMatrix=SILAC_XICs01_sepc_pep(IntervalScan(1):IntervalScan(2),5:8);
    LightMatrix_w=Matrix_Weight_func(LightMatrix);
    HeavyMatrix_w=Matrix_Weight_func(HeavyMatrix);
    HLRATIO_w(kkk)=sum(sum(HeavyMatrix_w))/sum(sum(LightMatrix_w));
    HLRATIO(kkk)=sum(sum(HeavyMatrix))/sum(sum(LightMatrix));
end
Unique_peptide_seq_string=unique(Peptide_seq);
for rrr=1:length(Unique_peptide_seq_string)
    J_Is=strcmp(Unique_peptide_seq_string{rrr},Peptide_seq);
    HLRATIO_w_uniquepepstr(rrr)=mean(HLRATIO_w(J_Is));
    HLRATIO_uniquepepstr(rrr)=mean(HLRATIO(J_Is));    
end
ID_noninf=~isinf(HLRATIO_w_uniquepepstr);
HLRATIO_w_uniquepepstrv1=HLRATIO_w_uniquepepstr(ID_noninf);
HLRATIO_uniquepepstrv1=HLRATIO_uniquepepstr(ID_noninf);
ID_nozero_w=find(HLRATIO_w_uniquepepstrv1~=0);
ID_nozero=find(HLRATIO_uniquepepstrv1~=0);
mean(log2(HLRATIO_w_uniquepepstrv1(ID_nozero_w)))
std(log2(HLRATIO_w_uniquepepstrv1(ID_nozero_w)))
mean(log2(HLRATIO_uniquepepstrv1(ID_nozero)))
std(log2(HLRATIO_uniquepepstrv1(ID_nozero)))
figure
hist(log2(HLRATIO_w_uniquepepstrv1(ID_nozero_w)),[-10:0.1:3])

HLRATIO=HLRATIO_w;
sum(isnan(HLRATIO))
sum(isinf(HLRATIO))
[Unique_Prot_id,III]=unique(Protein_hitnum);
Unique_Prot=Protein_seq(III);
id=1;id_up=1;
for i=1:length(Unique_Prot_id)
    ID_Same_prot=find(Protein_hitnum==Unique_Prot_id(i));
    Spec_Protein_unique_id=Protein_unique_id(ID_Same_prot);
    Spec_Peptide_seq_withmod=Peptide_seq_withmod(ID_Same_prot);
    Spec_Peptide_seq=Peptide_seq(ID_Same_prot);
    Spec_HLRATIO=HLRATIO(ID_Same_prot);
    
    ID_ninf=~isinf(Spec_HLRATIO);
    Spec_HLRATIOv1=Spec_HLRATIO(ID_ninf);
    Spec_Protein_unique_idv1=Spec_Protein_unique_id(ID_ninf);
    Spec_HLRATIOv2=Spec_HLRATIOv1(Spec_HLRATIOv1~=0);
    Spec_Protein_unique_idv2=Spec_Protein_unique_idv1(Spec_HLRATIOv1~=0);
    
    if ~isempty(Spec_HLRATIOv2)
        Prot_HLR(id)=(mean(unique(Spec_HLRATIOv2)));
        Prot_HLR_inter{id}=Spec_HLRATIOv2;
        Prot_HLR_inter_unique_pep{id}=Spec_HLRATIOv2(logical(Spec_Protein_unique_idv2));
%         Prot_HLR_inter_diff{id}=log2(Spec_HLRATIOv2)-mean(log2(unique(Spec_HLRATIOv2)));
        Prot_hitnum(id)=Unique_Prot_id(i);
        Prot_seq{id}=Unique_Prot{i};
        id=id+1;
        ID_uni=find(Spec_Protein_unique_idv2==1);
        if ~isempty(ID_uni)
            Prot_HLR_based_uniquepep(id_up)=mean(unique(Spec_HLRATIOv2(ID_uni)));
            Prot_hitnum_uniquepep(id)=Unique_Prot_id(i);
            Prot_seq_uniquepep{id}=Unique_Prot{i};
            id_up=id_up+1;
        end
    end
end
length(Prot_HLR)
figure
hist(Prot_HLR,[0:0.1:10])
ID_nonan=~isnan(Prot_HLR);
sum(ID_nonan)
mean(log2(Prot_HLR(ID_nonan)))
std(log2(Prot_HLR(ID_nonan)))
length(log2(Prot_HLR_based_uniquepep))
mean(log2(Prot_HLR_based_uniquepep))
std(log2(Prot_HLR_based_uniquepep))

Prot_HLR;
Prot_hitnum;
Prot_seq;
length(Protein_SC10_LHR)
length(Protein__SC10_seq_withLHR)
mean(log2(1./Protein_SC10_LHR))
std(log2(1./Protein_SC10_LHR))
Protein_SC10_hitnum;

num=1;num_uni=1;
for i=1:length(Prot_HLR_inter)
    if length(Prot_HLR_inter{i})>0
        Prot_mean(num)=mean(Prot_HLR_inter{i});
        Prot_median(num)=median(Prot_HLR_inter{i});
        Prot_mean_log2(num)=mean(log2(Prot_HLR_inter{i}));
        Prot_median_log2(num)=median(log2(Prot_HLR_inter{i}));
        Prot_std(num)=std(Prot_HLR_inter{i});
        num=num+1;
    end
    if length(Prot_HLR_inter_unique_pep{i})>0
        Prot_mean_uni_pep(num_uni)=mean(Prot_HLR_inter_unique_pep{i});
        Prot_median_uni_pep(num_uni)=median(Prot_HLR_inter_unique_pep{i});
        Prot_mean_log2_uni_pep(num_uni)=mean(log2(Prot_HLR_inter_unique_pep{i}));
        Prot_median_log2_uni_pep(num_uni)=median(log2(Prot_HLR_inter_unique_pep{i}));
        Prot_std_uni_pep(num_uni)=std(Prot_HLR_inter_unique_pep{i});
        num_uni=num_uni+1;
    end
end
ID_nonan=~isnan(Prot_mean);
mean(Prot_mean(ID_nonan))
mean(log2(Prot_mean(ID_nonan)))
mean(Prot_median)
mean(log2(Prot_median))
mean(Prot_mean_log2)
mean(Prot_median_log2)
std(log2(Prot_mean))
figure
hist(Prot_mean)

mean(Prot_mean_uni_pep)
mean(log2(Prot_mean_uni_pep))
mean(Prot_median_uni_pep)
mean(log2(Prot_median_uni_pep))
mean(Prot_mean_log2_uni_pep)
mean(Prot_median_log2_uni_pep)
std(log2(Prot_mean_uni_pep))


% [Com_hitnum,IA,IB]=intersect(Prot_hitnum,Protein_SC10_hitnum);
% n=9;
% Prot_hitnum(IA(n))
% Protein_SC10_hitnum(IB(n))
% Prot_seq(IA(n))
% Protein__SC10_seq_withLHR(IB(n))


Cal_HLR=[];
MAscort_HLR=[];
Cal_HLR_cmpMQ=[];        
Maxquant_HLR=[];
for i=1:length(Prot_seq)
    Spec_Protein_seq_Cal=Prot_seq{i};
    JUD_id=strcmp(Spec_Protein_seq_Cal,Protein__SC10_seq_withLHR);
    JUD_id_maxquant=strcmp(Spec_Protein_seq_Cal,Final_Maxquant_Razor_Protein_unique);
    if sum(JUD_id)>=1
        Cal_HLR=[Cal_HLR;Prot_HLR(i)];        
        MAscort_HLR=[MAscort_HLR;1./Protein_SC10_LHR(JUD_id)];
    end
    if sum(JUD_id_maxquant)>=1
        Cal_HLR_cmpMQ=[Cal_HLR_cmpMQ;Prot_HLR(i)];        
        Maxquant_HLR=[Maxquant_HLR;Final_Maxquant_Razor_Protein_unique_HLR(JUD_id_maxquant)];
    end   
    
end
ID_nonan=~isnan(Cal_HLR);
sum(ID_nonan)
length(log2(Cal_HLR(ID_nonan)))
mean(log2(Cal_HLR(ID_nonan)))
std(log2(Cal_HLR(ID_nonan)))
mean(log2(MAscort_HLR(ID_nonan)))
std(log2(MAscort_HLR(ID_nonan)))
figure
subplot(2,1,1);hist(log2(Cal_HLR(ID_nonan)),[-20:0.1:10]);grid on
subplot(2,1,2);hist(log2(MAscort_HLR(ID_nonan)),[-20:0.1:10]);grid on
ID_cal_good=find(log2(Cal_HLR(ID_nonan))>-2.5);
Mascort_ID_cal_good=find(log2(MAscort_HLR(ID_nonan))>-2.5);
Intersect_cal_Mascort=intersect(ID_cal_good,Mascort_ID_cal_good);
length(ID_cal_good)
length(Mascort_ID_cal_good)
length(Intersect_cal_Mascort)

ID_nonan=~isnan(Cal_HLR_cmpMQ);
sum(ID_nonan)
length(log2(Cal_HLR_cmpMQ(ID_nonan)))
mean((Cal_HLR_cmpMQ(ID_nonan)))
std(log2(Cal_HLR_cmpMQ(ID_nonan)))
mean(log2(Maxquant_HLR(ID_nonan)))
std(log2(Maxquant_HLR(ID_nonan)))
figure
subplot(2,1,1);hist(log2(Cal_HLR_cmpMQ(ID_nonan)),[-20:0.1:10]);grid on
subplot(2,1,2);hist(log2(Maxquant_HLR(ID_nonan)),[-20:0.1:10]);grid on
ID_cal_good=find(log2(Cal_HLR_cmpMQ(ID_nonan))>-2.5);
Maxquant_ID_cal_good=find(log2(Maxquant_HLR(ID_nonan))>-12.5);
Intersect_cal_Maxquant=intersect(ID_cal_good,Maxquant_ID_cal_good);
length(ID_cal_good)
length(Maxquant_ID_cal_good)
length(Intersect_cal_Maxquant)

sum(Protein_unique_id)
Unique_pepwithmod=unique(Peptide_seq_withmod);
Unique_pep=unique(Peptide_seq);
num=1;
for i=1:length(Unique_pep)
    Pep_u_seq=Unique_pep{i};
    ID_same=strcmp(Peptide_seq,Pep_u_seq);
    PEP_ch=Peptide_seq_withmod(ID_same);
    Prot_unique_ch=Protein_unique_id(ID_same);
    Prot_hitnum_ch=Protein_hitnum(ID_same);
    HLR_ch=HLRATIO(ID_same);
    PEP_ch_unique=unique(PEP_ch);
    for j=1:length(PEP_ch_unique)
        ID_same_ch=strcmp(PEP_ch,PEP_ch_unique{j});
        HLR_Final(num)=sum(HLR_ch(ID_same_ch))/length(HLR_ch(ID_same_ch));
        PEP_unique_Final{num}=PEP_ch_unique{j};
        Prot_unique_Final(num)=sum(Prot_unique_ch(ID_same_ch))/length(Prot_unique_ch(ID_same_ch));
        Prot_hitnum_ch_Final(num)=sum(Prot_hitnum_ch(ID_same_ch))/length(Prot_hitnum_ch(ID_same_ch));
        num=num+1;
    end    
end
ppp=2;
PEP_unique_Final{ppp}
Prot_unique_Final(ppp)
Prot_hitnum_ch_Final(ppp)
unique(Prot_hitnum_ch_Final)
num=1;
for k=min(Prot_hitnum_ch_Final):max(Prot_hitnum_ch_Final)
    ID_same=find(Prot_hitnum_ch_Final==k);
    if ~isempty(ID_same)
        HLR_spec=HLR_Final(ID_same);
        ID_ninf=~isinf(HLR_spec);
        if ~isempty(ID_ninf)
            HLR_specv1=HLR_spec(ID_ninf);
            ID_nnan=~isnan(HLR_specv1);
            if ~isempty(ID_nnan)
                HLR_Prot(num)=mean(log2(HLR_specv1(ID_nnan)));
                num=num+1;
            end
        end
    end
end
ID_ninf=~isinf(HLR_Prot);
mean((HLR_Prot(ID_ninf)))
std((HLR_Prot(ID_ninf)))
figure
hist(HLR_Prot(ID_ninf))
%%%%%%%%%%%%%%%%%%%%%%%%%

Num_protein_unique=0; ID_protein_unique=[];
for i=1:length(IntervalList01)
    if IntervalList01{i}.Labelling_efficiency_after7_combine(1)~=0 && Peptide_SC10_information.pep_isunique{i}==1       
        Num_protein_unique=Num_protein_unique+1;
        ID_protein_unique=[ID_protein_unique; i];
%         if IntervalList01{i}.Labelling_efficiency_after7_combine~=0 && ~isinf(IntervalList01{i}.Labelling_efficiency_after7_combine)    
%             num=num+1;
%             HLRatio=[HLRatio,IntervalList01{i}.Labelling_efficiency_after7_combine];
%         end    
    end
end
[U_pep,U_id]=unique(pep01(ID_protein_unique),'first');
sort(U_id);
ID_uniquepep_uniqueprot=sort(ID_protein_unique(U_id));
for i=1:length(ID_uniquepep_uniqueprot)
    row_id=ID_uniquepep_uniqueprot(i);
    
    
end

HLRatio=[];num=0;
for i=1:length(IntervalList01)    
    if length(IntervalList01{i}.Labelling_efficiency_after7_combine)==1        
        if IntervalList01{i}.Labelling_efficiency_after7_combine~=0 && ~isinf(IntervalList01{i}.Labelling_efficiency_after7_combine)    
            num=num+1;
            HLRatio=[HLRatio,IntervalList01{i}.Labelling_efficiency_after7_combine];
        end    
    end
end
figure
hist(log2(HLRatio),[-10:0.1:10])
grid on
std(log2(HLRatio))
mean(log2(HLRatio))
%%%%%%%%%%%%%%%%%%
ii=5569;
jj=4;

eval(['load ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_01_f', num2str(jj), '.mzXML.centroid0.peak.mat'])
Orbit_data01l1=peakl;
retentiont01_l1l2=retentiont;
id_level1_01=find(retentiont01_l1l2(:,2)==1);
retentiont01l1=retentiont01_l1l2(id_level1_01,1);
clear peakl retentiont

eval(['load ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_02_f', num2str(jj), '.mzXML.centroid0.peak.mat'])
Orbit_data02l1=peakl;
retentiont02_l1l2=retentiont;
id_level1_02=find(retentiont02_l1l2(:,2)==1);
retentiont02l1=retentiont02_l1l2(id_level1_02,1);
clear peakl retentiont

eval(['load ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_03_f', num2str(jj), '.mzXML.centroid0.peak.mat'])
Orbit_data03l1=peakl;
retentiont03_l1l2=retentiont;
id_level1_03=find(retentiont03_l1l2(:,2)==1);
retentiont03l1=retentiont03_l1l2(id_level1_03,1);
clear peakl retentiont
    
clear data01l1 data02l1 data03l1
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
load AMT_database
%%%%%%%%%%%
Total_time_diff=[];
for i=1:length(AMT_database)
    Time_VEC=[];
   for j=1:length(AMT_database(i).ortherinformation)
        Time_VEC=[Time_VEC;AMT_database(i).ortherinformation{j}.normal_rt];
   end
   Total_time_diff=[Total_time_diff;Time_VEC-mean(Time_VEC)];
end
var(Total_time_diff);
Std_AMT_timeshift=std(Total_time_diff);
save D:\Program\QTOF_replicate_identification\MATfile_v2\Std_AMT_timeshift Std_AMT_timeshift
% %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Generate AMT data base on fraction 5569 f4
% load Align_Out_matrix_5569_f4 %%%% based on long txt file
load XUEPO_Orbi_5569_f4_info %%%% based on xuepo xlcs file
load TOTAL_XICs_Orbitrap_5569_f4%%% Orbi XICs

for i=20:25
    colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
    S_start=str2num(matrix{i,6});
    S_end=str2num(matrix{i,7});
    j=posi01(str2num(matrix{i,9}));
    figure
    for k=1:6
        plot(XICs_Orbi01{k}(:,j),colorarray(k))
        hold on
    end   
    stem(S_start,10^7,'ro')
    stem(S_end,10^7,'ko')
end

Orbit_Posi=zeros(size(matrix,1)-1,3);
for i=2:size(matrix,1)
    Orbit_Posi(i-1,1)=str2num(matrix{i,9});
    Orbit_Posi(i-1,2)=str2num(matrix{i,17});
    Orbit_Posi(i-1,3)=str2num(matrix{i,25});
    Orbit_time(i-1,1)=str2num(matrix{i,4});
    Orbit_time(i-1,2)=str2num(matrix{i,6});
    Orbit_time(i-1,3)=str2num(matrix{i,7});
    Orbit_time(i-1,4)=str2num(matrix{i,12});
    Orbit_time(i-1,5)=str2num(matrix{i,14});
    Orbit_time(i-1,6)=str2num(matrix{i,15});
    Orbit_time(i-1,7)=str2num(matrix{i,20});
    Orbit_time(i-1,8)=str2num(matrix{i,22});
    Orbit_time(i-1,9)=str2num(matrix{i,23});
    Orbit_cs(i-1,1)=str2num(matrix{i,2});
    Orbit_cs(i-1,2)=str2num(matrix{i,10});
    Orbit_cs(i-1,3)=str2num(matrix{i,18});
end

for k=1:size(Orbit_Posi,1)
    Jud_vect=Orbit_Posi(k,:)~=0;
    Jud_num=sum(Jud_vect.*[1 2 4]);
    switch Jud_num
        case {1,3,5,7} 
            for i=1:6
                Orbit_TotalXICs01{i}(:,k)=XICs_Orbi01{i}(:,posi01(Orbit_Posi(k,1)));
                Orbit_TotalXICs02{i}(:,k)=XICs_Orbi0102{i}(:,posi01(Orbit_Posi(k,1)));
                Orbit_TotalXICs03{i}(:,k)=XICs_Orbi0103{i}(:,posi01(Orbit_Posi(k,1)));
            end
        case {2,6}
            for i=1:6
                Orbit_TotalXICs01{i}(:,k)=XICs_Orbi0201{i}(:,posi02(Orbit_Posi(k,2)));
                Orbit_TotalXICs02{i}(:,k)=XICs_Orbi02{i}(:,posi02(Orbit_Posi(k,2)));
                Orbit_TotalXICs03{i}(:,k)=XICs_Orbi0203{i}(:,posi02(Orbit_Posi(k,2)));
            end
        case {4}
            for i=1:6
                Orbit_TotalXICs01{i}(:,k)=XICs_Orbi0301{i}(:,posi03(Orbit_Posi(k,3)));
                Orbit_TotalXICs02{i}(:,k)=XICs_Orbi0302{i}(:,posi03(Orbit_Posi(k,3)));
                Orbit_TotalXICs03{i}(:,k)=XICs_Orbi03{i}(:,posi03(Orbit_Posi(k,3)));
            end
    end
end
  clear   XICs_Orbi01 XICs_Orbi02 XICs_Orbi03 XICs_Orbi0102 XICs_Orbi0103 XICs_Orbi0201 XICs_Orbi0203 XICs_Orbi0301 XICs_Orbi0302
for i=2:size(matrix,1)
    
    AMT_database_5569_f4(i-1).peptide=matrix{i,1};
    
    CS_vector=[str2num(matrix{i,2}),str2num(matrix{i,10}),str2num(matrix{i,18})];
    Jud_vec=find(CS_vector~=0);       
    AMT_database_5569_f4(i-1).charge=CS_vector(Jud_vec(1));
    
    Time_vector=[str2num(matrix{i,4}),str2num(matrix{i,12}),str2num(matrix{i,20})];
    Jud_vec=find(Time_vector~=0);       
    AMT_database_5569_f4(i-1).normal_rt=sum(Time_vector(Jud_vec))/length(Jud_vec);
    
    [peptidenew, massdiffList, isHeavy]=modprocess({AMT_database_5569_f4(i-1).peptide});
    peptidenew=peptidenew{1};
    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(AMT_database_5569_f4(i-1).peptide));
    AMT_database_5569_f4(i-1).mass=weight;
    csvalue=AMT_database_5569_f4(i-1).charge;
    mzvalue=(weight+csvalue*1.0073)/csvalue;
    AMT_database_5569_f4(i-1).mono_mz=mzvalue;
    AMT_database_5569_f4(i-1).iso=isotopepattern;
end

%%%%%%%%% calculate labeling effeciency
for i=1:length(AMT_database_5569_f4)
    isoList=AMT_database_5569_f4(i).iso;
    MS2detecttime_vector=[Orbit_time(i,1),Orbit_time(i,4),Orbit_time(i,7)];
    MS2detectID_vector=find(MS2detecttime_vector~=0);
    MS2detectID=MS2detectID_vector(1);
    
    Scan_detect=Orbit_time(i,(MS2detectID-1)*3+2:(MS2detectID-1)*3+3);
    for j=1:6
        EV=['MS2_Intervaldata(:,j)=Orbit_TotalXICs0',num2str(MS2detectID),'{j}(Scan_detect(1):Scan_detect(2),i);'];
        eval(EV);
    end
    
    minf=0.0;maxf=1;rangerate=0.6;
    intensity=sum(MS2_Intervaldata,1);
    est=O18rateLinear(isoList(1:6),intensity,minf,maxf,rangerate,1);
    Orbit_Labelling_efficiency(i)=est(3);
    Orbit_Spec_Pep_O18rate(i)=est(2);
    clear MS2_Intervaldata   
    
end
%%%%%%%%%

%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% load TOF raw data
filename=['ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_01_f4.mzXML.centroid0.peak.mat'];
load(filename);
TOF_peakl01=peakl;
%%%%%%%%%%% cut off the repeat part
ID_turn01=find(retentiont(1:end-1)-retentiont(2:end)>0);
TOF_retentiont01l1_orignal=retentiont;
TOF_retentiont01l1=retentiont(1:ID_turn01(1));
%%%%%%%%%%%
clear peakl filename retentiont;

filename=['ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_02_f4.mzXML.centroid0.peak.mat'];
load(filename);
TOF_peakl02=peakl;
%%%%%%%%%%% cut off the repeat part
ID_turn02=find(retentiont(1:end-1)-retentiont(2:end)>0);
TOF_retentiont02l1_orignal=retentiont;
TOF_retentiont02l1=retentiont(1:ID_turn02(1));
%%%%%%%%%%%
clear peakl filename retentiont;

filename=['ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_03_f4.mzXML.centroid0.peak.mat'];
load(filename);
TOF_peakl03=peakl;
%%%%%%%%%%% cut off the repeat part
ID_turn03=find(retentiont(1:end-1)-retentiont(2:end)>0);
TOF_retentiont03l1_orignal=retentiont;
TOF_retentiont03l1=retentiont(1:ID_turn03(1));
%%%%%%%%%%%
clear peakl filename retentiont;
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
Record_AMT_same_TOFv1=[1:length(AMT_database_5569_f4)]';
TOFAMT_LIST=[];
for i=1:length(AMT_database_5569_f4)
    ID_AMT=find(Record_AMT_same_TOFv1(:,1)==i);
    if ~isempty(ID_AMT)
        TOFAMT_LIST=[TOFAMT_LIST;Record_AMT_same_TOFv1(ID_AMT(1),:)];  
    end
end

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% get XICs with long's function
for i=1:size(TOFAMT_LIST,1)
    AMT_database_5569_f4_inLongTOFpep(i)=AMT_database_5569_f4(TOFAMT_LIST(i,1));
end
Orbit_AMTms2time=Orbit_time(TOFAMT_LIST(:,1),:);
Orbit_AMTcs=Orbit_cs(TOFAMT_LIST(:,1),:);
for i=1:6
    Orbit_AMTXICs01{i}=Orbit_TotalXICs01{i}(:,TOFAMT_LIST(:,1));
    Orbit_AMTXICs02{i}=Orbit_TotalXICs02{i}(:,TOFAMT_LIST(:,1));
    Orbit_AMTXICs03{i}=Orbit_TotalXICs03{i}(:,TOFAMT_LIST(:,1));
end

TOF_totalmzList=[];
for i=1:length(AMT_database_5569_f4_inLongTOFpep)
    [peptidenew, massdiffList, isHeavy]=modprocess({AMT_database_5569_f4_inLongTOFpep(i).peptide});
    peptidenew=peptidenew{1};
    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(AMT_database_5569_f4_inLongTOFpep(i).peptide));
    csvalue=AMT_database_5569_f4_inLongTOFpep(i).charge;
    mzvalue=(weight+csvalue*1.0073)/csvalue;
    TOF_totalmzList=[TOF_totalmzList; mzvalue mzvalue+1.00335/csvalue mzvalue+2.00547/csvalue mzvalue+3.00882/csvalue mzvalue+4.00849/csvalue mzvalue+5.01184/csvalue];
end

tolerance=80;
for i=1:6
    TOF_XICs01_Tolerance80{i}=getXIC_LC_new(peakl01,TOF_totalmzList(:,i),tolerance);
    TOF_XICs02_Tolerance80{i}=getXIC_LC_new(peakl02,TOF_totalmzList(:,i),tolerance);
    TOF_XICs03_Tolerance80{i}=getXIC_LC_new(peakl03,TOF_totalmzList(:,i),tolerance);
end
clear tolerance
save D:\Program\QTOF_replicate_identification\MATfile_v2\TOF_XICS_tolerance80 TOF_XICs01_Tolerance80 TOF_XICs02_Tolerance80 TOF_XICs03_Tolerance80
load D:\Program\QTOF_replicate_identification\MATfile_v2\TOF_XICS_tolerance80 TOF_XICs01_Tolerance80 TOF_XICs02_Tolerance80 TOF_XICs03_Tolerance80
for i=1:6
    TOF_XICs01{i}=TOF_XICs01_Tolerance80{i}(1:length(TOF_retentiont01l1),:);
    TOF_XICs02{i}=TOF_XICs02_Tolerance80{i}(1:length(TOF_retentiont02l1),:);
    TOF_XICs03{i}=TOF_XICs03_Tolerance80{i}(1:length(TOF_retentiont03l1),:);
end
clear TOF_XICs01_Tolerance80 TOF_XICs02_Tolerance80 TOF_XICs03_Tolerance80
% tolerance=20;
% for i=1:6
%     TOF_XICs01_Tolerance20{i}=getXIC_LC_new(peakl01,TOF_totalmzList(:,i),tolerance);
%     TOF_XICs02_Tolerance20{i}=getXIC_LC_new(peakl02,TOF_totalmzList(:,i),tolerance);
%     TOF_XICs03_Tolerance20{i}=getXIC_LC_new(peakl03,TOF_totalmzList(:,i),tolerance);
% end
% clear tolerance
% save D:\Program\QTOF_replicate_identification\MATfile_v2\TOF_XICS_tolerance20 TOF_XICs01_Tolerance20 TOF_XICs02_Tolerance20 TOF_XICs03_Tolerance20
% 
% for i=1:6
%     TOF_XICs01_20ppm{i}=TOF_XICs01_Tolerance20{i}(1:length(TOF_retentiont01l1),:);
%     TOF_XICs02_20ppm{i}=TOF_XICs02_Tolerance20{i}(1:length(TOF_retentiont02l1),:);
%     TOF_XICs03_20ppm{i}=TOF_XICs03_Tolerance20{i}(1:length(TOF_retentiont03l1),:);
% end
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% interval detection
for i=1:length(AMT_database_5569_f4_inLongTOFpep)
    Pepcommon0102v1{i}=AMT_database_5569_f4_inLongTOFpep(i).peptide;
    TOF_ms2information01(i)=AMT_database_5569_f4_inLongTOFpep(i).normal_rt;
    TOF_ms2information02=TOF_ms2information01;
    TOF_ms2information03=TOF_ms2information01;
    iso(i,:)=AMT_database_5569_f4_inLongTOFpep(i).iso;
    TOF_mass(i)=AMT_database_5569_f4_inLongTOFpep(i).mass;
    TOF_charge(i)=AMT_database_5569_f4_inLongTOFpep(i).charge;
end

Iso12KL_th=Generate_iso12KL_th(iso);

pep01=Pepcommon0102v1;
posi01v1=[];
posi01v2=[];
posi01=[];
for i=1:length(pep01)
    TOF_XICs01_sepc_pep=zeros(size(TOF_XICs01{1},1),6);
    for j=1:6
        TOF_XICs01_sepc_pep(:,j)=TOF_XICs01{j}(:,i);        
    end
    [IntervalList,Jud_detect_good]=TOF_Verification_intervaldetection_Spec_Pep(pep01{i},iso(i,:),TOF_ms2information01(i),TOF_XICs01_sepc_pep,TOF_retentiont01l1);
    IntervalList01{i}=IntervalList;
    if Jud_detect_good==1
        posi01v1=[posi01v1;i];
        if IntervalList.intervallist_after7_combine(1,1)~=0 && IntervalList.intervallist_after7_combine(1,2)~=0
            posi01v2=[posi01v2;i];
%             TOF_retentiont01l1v1=TOF_retentiont01l1';
            T_m=TOF_retentiont01l1(IntervalList.intervallist_after7_combine)-IntervalList.ms2time;
            T_v=T_m(:,1).*T_m(:,2);
            Id_ms2interval=find(T_v<=0);
            if ~isempty(Id_ms2interval)
                posi01=[posi01;i];
            end            
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
%             TOF_retentiont02l1v1=TOF_retentiont02l1';
            T_m=TOF_retentiont02l1(IntervalList.intervallist_after7_combine)-IntervalList.ms2time;
            T_v=T_m(:,1).*T_m(:,2);
            Id_ms2interval=find(T_v<=0);
            if ~isempty(Id_ms2interval)
                posi02=[posi02;i];
            end            
        end
    end    
end

save D:\Program\QTOF_replicate_identification\MATfile_v2\intervallist IntervalList01 IntervalList02

length(posi01v1)
length(posi01v2)
length(posi01)
length(posi02v1)
length(posi02v2)
length(posi02)

save D:\Program\QTOF_replicate_identification\MATfile_v2\IDs posi01v1 posi01v2 posi01 posi02v1 posi02v2 posi02

%%%%%%%%%%%%%%% check interval detection
%%%%%% check iso(1:2) to Orbit ms2 interval
ID_ms2detect_Orbit01=find(Orbit_AMTms2time(:,1)~=0);
for k=1:length(ID_ms2detect_Orbit01)
    i=ID_ms2detect_Orbit01(k);
    for j=1:6
        Orbit_inteval(:,j)=Orbit_TotalXICs01{j}(Orbit_AMTms2time(i,2):Orbit_AMTms2time(i,3),i);
    end
    Orbit_iso12_pattern=sum(Orbit_inteval(:,1:2),1)./sum(sum(Orbit_inteval(:,1:2)));
    log_iso12_KL_Value(k)=KL_calculate(Orbit_iso12_pattern,iso(i,1:2));    
    clear Orbit_inteval
end
figure
hist(log_iso12_KL_Value,100)

Sorted_log_sio12_KL_value=sort(log_iso12_KL_Value);
Iso12KL_th=Sorted_log_sio12_KL_value(round(length(Sorted_log_sio12_KL_value)*0.95));
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
Times_noise_std=3;%6;
instrument_h_int=10^7;
de_range=10^5;%10^4;
colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
for i=36:40%%length(posi01v2)
    ID_into=posi01v2(i);
    Labeling_eff_diff=abs(IntervalList01{ID_into}.Orbit_Labelling_efficiency_after7_combine-Orbit_Labelling_efficiency(ID_into));
%     for j=1:6
%         th(j)=TOF_getThreshold_range(TOF_XICs01{j}(:,ID_into),instrument_h_int, de_range,Times_noise_std);
%         intervalRawData(:,j)=TOF_XICs01{j}(IntervalList01{ID_into}.intervallist_after7_combine(2,1):IntervalList01{ID_into}.intervallist_after7_combine(2,2),ID_into);
%     end
%     Judge=long_criteria(iso(ID_into,:),intervalRawData,th);
    
    figure
    subplot(2,2,1)
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
    for kkk=1:length(IntervalList01{ID_into}.intervallist_after7_combine(:,1))
        interval_mid_time=(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(kkk,1))+TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(kkk,2)))/2;
        text(interval_mid_time,2*Height,num2str(Labeling_eff_diff(kkk)));
    end
    grid on
    subplot(2,2,2)
    for j=1:6
        plot(retentiont01l1,Orbit_TotalXICs01{j}(:,ID_into),colorarray(j))
        hold on
    end
    stem(retentiont01l1(Orbit_AMTms2time(ID_into,2)),10^7,'ro')
    stem(retentiont01l1(Orbit_AMTms2time(ID_into,3)),10^7,'ko')
    grid on    
    
    subplot(2,2,3)
    for j=1:6
        plot(TOF_retentiont01l1,TOF_XICs01_20ppm{j}(:,ID_into),colorarray(j))
        hold on
    end
    grid on
    
%     MZ_interval(1)=TOF_totalmzList(ID_into,1)-0.2;
%     MZ_interval(2)=TOF_totalmzList(ID_into,6)+0.2;
%     Interval_start=IntervalList01{ID_into}.intervallist_after7_combine(1,1);
%     Interval_end=IntervalList01{ID_into}.intervallist_after7_combine(1,2);
%     [Int_Matrix_mz_time,mz_unique,Pep_mz,Interval_XIC]=Generate_Pep_Mz_Time_matrix(TOF_peakl01,MZ_interval,Interval_start,Interval_end);
%     x=mz_unique*ones(1,size(Int_Matrix_mz_time,2));
%     y=ones(size(Int_Matrix_mz_time,1),1)*[Interval_start:Interval_end];
%     figure
%     surf(x,y,Int_Matrix_mz_time)    
%     figure
%     plot(mz_unique,sum(Int_Matrix_mz_time,2))
%     hold on
%     grid on
%     stem(TOF_totalmzList(ID_into,:),30000*ones(1,6),'ro')
    
end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%

[Com_ID_0102,ID_0102_01,ID_0102_02]=intersect(posi01v2,posi02v2);
Orbit_AMTms2time_com=Orbit_AMTms2time(posi01v2(ID_0102_01),:);
Orbit_AMTcs_com=Orbit_AMTcs(posi01v2(ID_0102_01),:);

for i=1:6
    Orbit_XICs_com01{i}=Orbit_AMTXICs01{i}(:,posi01v2(ID_0102_01));
    Orbit_XICs_com02{i}=Orbit_AMTXICs02{i}(:,posi02v2(ID_0102_02));
    TOF_XICs_com01{i}=TOF_XICs01{i}(:,posi01v2(ID_0102_01));
    TOF_XICs_com02{i}=TOF_XICs02{i}(:,posi02v2(ID_0102_02));
end
IntervalList_com01=IntervalList01(posi01v2(ID_0102_01));
IntervalList_com02=IntervalList02(posi02v2(ID_0102_02));

Pep_com=Pepcommon0102v1(Com_ID_0102);
TOF_ms2information_com=TOF_ms2information01(Com_ID_0102);
iso_com=iso(Com_ID_0102,:);
TOF_mass_com=TOF_mass(Com_ID_0102);
TOF_charge_com=TOF_charge(Com_ID_0102);
Orbit_Labelling_efficiency_com=Orbit_Labelling_efficiency(Com_ID_0102);

save D:\Program\QTOF_replicate_identification\MATfile_v2\TOF_Total_information_cent0 Pep_com TOF_ms2information_com iso_com TOF_XICs_com01 TOF_XICs_com02 IntervalList_com01 IntervalList_com02
save D:\Program\QTOF_replicate_identification\MATfile_v2\Orbit_Total_information Pep_com Orbit_AMTms2time_com iso_com Orbit_XICs_com01 Orbit_XICs_com02 Orbit_Labelling_efficiency_com
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Generate all the AT AR AKL and R T LE scores(not probability) for
%%%%%%%%%%%%%%%%%% only TOF data
Orbit_AMTms2time_com;
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
    
    %%%%%%%%%%%%%%%%%%%% TOF01 to Orbit (R T LE)
    %%%%%%%%%%%%%%%%%%%%
    Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
    Interval_matrix02=Orbit_AMTms2time_com(i,2:3);
    if sum(Interval_matrix02)~=0
        TOFscorev1{i}.LE_TOF01Orbit01=IntervalList_com01{i}.Labelling_efficiency_after7_combine'-Orbit_Labelling_efficiency_com(i);        
        TOFscorev1{i}.LE_TOF01=IntervalList_com01{i}.Labelling_efficiency_after7_combine';
        for k=1:size(Interval_matrix01,1)
            for j=1:6
                XICs01{k}(:,j)=TOF_XICs_com01{j}(Interval_matrix01(k,1):Interval_matrix01(k,2),i);
            end 
        end
        for k=1:size(Interval_matrix02,1)
            for j=1:6
                XICs02{k}(:,j)=Orbit_XICs_com01{j}(Interval_matrix02(k,1):Interval_matrix02(k,2),i);
            end 
        end
%         for j=1:6
%             XICs01(:,j)=TOF_XICs_com01{j}(:,i);
%             XICs02(:,j)=Orbit_XICs_com01{j}(:,i);
%         end
        %%%%%% intervals01 and intervals02 are matrix;
        %%%%%% XICs01 and XICs02 are matrixs (n*6) 6 colums are for different
        %%%%%% iso mz
        [TOFscorev1{i}.LOGKL_TOF01Orbit01,TOFscorev1{i}.R_TOF01Orbit01,TOFscorev1{i}.T_TOF01Orbit01,TOFscorev1{i}.Normal_T_TOF01Orbit01,TOFscorev1{i}.T_TOF01,TOFscorev1{i}.Normal_T_TOF01,TOFscorev1{i}.T_Orbit01,TOFscorev1{i}.Normal_T_Orbit01]=GenerateRTKL(Interval_matrix01,Interval_matrix02,XICs01,XICs02,TOF_retentiont01l1,retentiont01l1,iso_spec_pep);
        
        clear XICs01 XICs02
    else
       TOFscorev1{i}.LOGKL_TOF01Orbit01=0;
       TOFscorev1{i}.R_TOF01Orbit01=0;
       TOFscorev1{i}.T_TOF01Orbit01=100000;
       TOFscorev1{i}.Normal_T_TOF01Orbit01=1;
       TOFscorev1{i}.T_TOF01=0;
       TOFscorev1{i}.Normal_T_TOF01=1;
       TOFscorev1{i}.T_Orbit01=0;
       TOFscorev1{i}.Normal_T_Orbit01=1;
       TOFscorev1{i}.LE_TOF01Orbit01=1;        
       TOFscorev1{i}.LE_TOF01=-1;
    end
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% TOF02 to Orbit
    %%%%%%%%%%%%%%%%%%%%
    Interval_matrix01=IntervalList_com02{i}.intervallist_after7_combine;
    Interval_matrix02=Orbit_AMTms2time_com(i,5:6);
    if sum(Interval_matrix02)~=0
        TOFscorev1{i}.LE_TOF02Orbit02=IntervalList_com02{i}.Labelling_efficiency_after7_combine'-Orbit_Labelling_efficiency_com(i);        
        TOFscorev1{i}.LE_TOF02=IntervalList_com02{i}.Labelling_efficiency_after7_combine';
        for k=1:size(Interval_matrix01,1)
            for j=1:6
                XICs01{k}(:,j)=TOF_XICs_com02{j}(Interval_matrix01(k,1):Interval_matrix01(k,2),i);
            end 
        end
        for k=1:size(Interval_matrix02,1)
            for j=1:6
                XICs02{k}(:,j)=Orbit_XICs_com02{j}(Interval_matrix02(k,1):Interval_matrix02(k,2),i);
            end 
        end
        
%         for j=1:6
%             XICs01(:,j)=TOF_XICs_com02{j}(:,i);
%             XICs02(:,j)=Orbit_XICs_com02{j}(:,i);
%         end
        %%%%%% intervals01 and intervals02 are matrix;
        %%%%%% XICs01 and XICs02 are matrixs (n*6) 6 colums are for different
        %%%%%% iso mz
        [TOFscorev1{i}.LOGKL_TOF02Orbit02,TOFscorev1{i}.R_TOF02Orbit02,TOFscorev1{i}.T_TOF02Orbit02,TOFscorev1{i}.Normal_T_TOF02Orbit02,TOFscorev1{i}.T_TOF02,TOFscorev1{i}.Normal_T_TOF02,TOFscorev1{i}.T_Orbit02,TOFscorev1{i}.Normal_T_Orbit02]=GenerateRTKL(Interval_matrix01,Interval_matrix02,XICs01,XICs02,TOF_retentiont02l1,retentiont02l1,iso_spec_pep);
        clear XICs01 XICs02
    else
       TOFscorev1{i}.LOGKL_TOF02Orbit02=0;
       TOFscorev1{i}.R_TOF02Orbit02=0;
       TOFscorev1{i}.T_TOF02Orbit02=100000;
       TOFscorev1{i}.Normal_T_TOF02Orbit02=1;
       TOFscorev1{i}.T_TOF02=0;
       TOFscorev1{i}.Normal_T_TOF02=1;
       TOFscorev1{i}.T_Orbit02=0;
       TOFscorev1{i}.Normal_T_Orbit02=1;
       TOFscorev1{i}.LE_TOF02Orbit02=1;
       TOFscorev1{i}.LE_TOF02=-1;
    end
    %%%%%%%%%%%%%%%%%%%%

end

save D:\Program\QTOF_replicate_identification\MATfile_v2\TOFscorev1 TOFscorev1
load D:\Program\QTOF_replicate_identification\MATfile_v2\TOFscorev1
TOFscore=TOFscorev1;
save D:\Program\QTOF_replicate_identification\MATfile_v2\TOFscore TOFscore

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% generate AR AT AKL model based on the training data
%%%%%%%%%%%%%%%%%% selected by T model and R model from TOF and Orbit
Training_Row_Column_Num=[];
Timepair_corr=[];
Timepair_corr_normal=[];
Timepair_noncorr=[];
Timepair_noncorr_normal=[];
AT_corr=[];AT_non_corr=[];
AR_corr=[];AR_non_corr=[];
AKL0102_corr=[];AKL0102_non_corr=[];
Training_id=[];
num1=0; num2=0; num3=0; num4=0;
for i=1:length(TOFscore)%length(Training_id_R_T)
%     i=Training_id_R_T(ii);
    Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
    Interval_matrix02=IntervalList_com02{i}.intervallist_after7_combine;
         
    Normalized_ms2time01_TOF=TOF_ms2information_com(i)/max(retentiont01l1);
    Normalized_ms2time02_TOF=TOF_ms2information_com(i)/max(retentiont02l1);
    
    R_AT=find(TOFscore{i}.Normal_T01(:,1)<=Normalized_ms2time01_TOF+2*Std_AMT_timeshift & TOFscore{i}.Normal_T01(:,1)>=Normalized_ms2time01_TOF-0*Std_AMT_timeshift);
    C_AT=find(TOFscore{i}.Normal_T02(1,:)<=Normalized_ms2time02_TOF+2*Std_AMT_timeshift & TOFscore{i}.Normal_T02(1,:)>=Normalized_ms2time02_TOF-0*Std_AMT_timeshift);
    %%%%% 1. Find any intervals in TOF that is in Orbit ms2 interval
    %%%%% 2. LOGAKL <= -2.5
    %%%%% 3. AR >= 0.85
%     if ~isempty(R_AT) && ~isempty(C_AT)
    %%%%% 1. Orbit ms2 interval only contains one interval in TOF
    %%%%% 2. LOGAKL <= -2.5
    %%%%% 3. AR >= 0.85
    if length(R_AT)==1 && length(C_AT)==1
        
        num1=num1+1;

        V_LOGAKL=TOFscore{i}.LOGAKL(R_AT,C_AT);
        V_AR=TOFscore{i}.AR(R_AT,C_AT);
        
        [AKL_R_good_id,AKL_C_good_id,AKL_V_good_id]=find(V_LOGAKL<=-2.5);
        [AR_R_good_id,AR_C_good_id,AR_V_good_id]=find(V_AR>=0.80);
        
        AKL_good_id_vector=[AKL_R_good_id,AKL_C_good_id];
        AR_good_id_vector=[AR_R_good_id,AR_C_good_id];
        
        Good_R_id=[];
        Good_C_id=[];
        if ~isempty(AKL_good_id_vector) && ~isempty(AR_good_id_vector)
            
            num2=num2+1;
            
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
                
                num3=num3+1;
                
                R_corr=R_AT(Good_R_id);
                C_corr=C_AT(Good_C_id);
                
                if length(R_corr)==1 && length(C_corr)==1

                    num4=num4+1;
                    
                    Training_id=[Training_id; i];

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

                    V_AR_corr=TOFscore{i}.AR(R_corr,C_corr);
                    Vector_AR_noncorr_Row=TOFscore{i}.AR(R_corr,:);
                    Vector_AR_noncorr_Col=TOFscore{i}.AR(:,C_corr);
                    Vector_AR_noncorr_Row(C_corr)=[];
                    Vector_AR_noncorr_Col(R_corr)=[];
                    AR_corr=[AR_corr;V_AR_corr];
                    AR_non_corr=[AR_non_corr;Vector_AR_noncorr_Row';Vector_AR_noncorr_Col];

                    V_AKL_corr=TOFscore{i}.LOGAKL(R_corr,C_corr);
                    Vector_AKL_noncorr_Row=TOFscore{i}.LOGAKL(R_corr,:);
                    Vector_AKL_noncorr_Col=TOFscore{i}.LOGAKL(:,C_corr);
                    Vector_AKL_noncorr_Row(C_corr)=[];
                    Vector_AKL_noncorr_Col(R_corr)=[];
                    AKL0102_corr=[AKL0102_corr;V_AKL_corr];
                    AKL0102_non_corr=[AKL0102_non_corr;Vector_AKL_noncorr_Row';Vector_AKL_noncorr_Col];        
                    
                    Training_Row_Column_Num=[Training_Row_Column_Num; R_corr C_corr];
                end
            end
        end
    end    
end
length(Training_id)
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% check training data
%%%%%%%%%%%%%%%
clear Interval_data01
clear log_KL_Value_corr
clear Interval_data01 Interval_data02
clear log_KL_Value_noncorr01
num=1;
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

num=1;
num_o=1;
ID_DEtect=[];
for i=1:length(Training_id)
    ID_into=Com_ID_0102(Training_id(i)); 
    Training_Pep{i}=Pep_com{Training_id(i)};
    if Orbit_AMTms2time(ID_into,2)~=0 && Orbit_AMTms2time(ID_into,3)~=0 && Orbit_AMTms2time(ID_into,5)~=0 && Orbit_AMTms2time(ID_into,6)~=0
        for j=1:6
            Orbit_training_data01.intervalsdata(:,j)=Orbit_TotalXICs01{j}(Orbit_AMTms2time(ID_into,2):Orbit_AMTms2time(ID_into,3),ID_into);
        end
        for j=1:6
            Orbit_training_data02.intervalsdata(:,j)=Orbit_TotalXICs02{j}(Orbit_AMTms2time(ID_into,5):Orbit_AMTms2time(ID_into,6),ID_into);
        end

        Orbit_training_data01.iso=iso_com(Training_id(i),:);
        Orbit_training_data02.iso=iso_com(Training_id(i),:);
        [Orbit_finalO18rate01(num_o,:),Orbit_finalO18f01(num_o,:)]=OrbitrapProduceO18RatesV1(Orbit_training_data01);
        [Orbit_finalO18rate02(num_o,:),Orbit_finalO18f02(num_o,:)]=OrbitrapProduceO18RatesV1(Orbit_training_data02);
        num_o=num_o+1;
        ID_DEtect=[ID_DEtect;i];
        clear Orbit_training_data01 Orbit_training_data02
%     end

        data01_intervalid=Training_Row_Column_Num(i,1);
        data02_intervalid=Training_Row_Column_Num(i,2);
        TOF_tr_scan_start01=IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,1);
        TOF_tr_scan_end01=IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,2);
        TOF_tr_scan_start02=IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,1);
        TOF_tr_scan_end02=IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,2);
        for j=1:6
            TOF_training_data01.intervalsdata(:,j)=TOF_XICs01{j}(TOF_tr_scan_start01:TOF_tr_scan_end01,ID_into);
        end   
        for j=1:6
            TOF_training_data02.intervalsdata(:,j)=TOF_XICs02{j}(TOF_tr_scan_start02:TOF_tr_scan_end02,ID_into);
        end
        TOF_training_data01.iso=iso_com(Training_id(i),:);
        TOF_training_data02.iso=iso_com(Training_id(i),:);

        [TOF_finalO18rate01(num,:),TOF_finalO18f01(num,:)]=OrbitrapProduceO18RatesV1(TOF_training_data01); 
        [TOF_finalO18rate02(num,:),TOF_finalO18f02(num,:)]=OrbitrapProduceO18RatesV1(TOF_training_data02);
        num=num+1;
        clear TOF_training_data01 TOF_training_data02
    end
        %%%%%%%%%%%%%%%%%%%
end
save Training_Pep Training_Pep ID_DEtect
std(log2(TOF_finalO18rate01(:,3)))
std(log2(TOF_finalO18rate02(:,3)))
std(log2(TOF_finalO18rate01(:,3))-log2(TOF_finalO18rate02(:,3)))
mean(log2(TOF_finalO18rate01(:,3))-log2(TOF_finalO18rate02(:,3)))
std(log2(Orbit_finalO18rate01(:,3)))
std(log2(Orbit_finalO18rate02(:,3)))
std(log2(Orbit_finalO18rate01(:,3))-log2(Orbit_finalO18rate02(:,3)))
mean(log2(Orbit_finalO18rate01(:,3))-log2(Orbit_finalO18rate02(:,3)))

std(log2(TOF_finalO18rate01(ID2,3))-log2(TOF_finalO18rate02(ID2,3)))
std(log2(Orbit_finalO18rate01(ID2,3))-log2(Orbit_finalO18rate02(ID2,3)))

std(log2(Orbit_finalO18rate01(:,3))-log2(TOF_finalO18rate01(:,3)))
mean(log2(Orbit_finalO18rate01(:,3))-log2(TOF_finalO18rate01(:,3)))
%%%%%%%%%%%% random range the id of the Orbit
% L=length(Orbit_finalO18rate01(:,3))
% ID_picked=randsample(length(Orbit_finalO18rate01(:,3)),L,'false');
% std(log2(Orbit_finalO18rate01(:,3))-log2(Orbit_finalO18rate02(ID_picked,3)))
% mean(log2(Orbit_finalO18rate01(:,3))-log2(Orbit_finalO18rate02(ID_picked,3)))
num=1;
for L=1:20:200;
    ID_picked=randsample(length(Orbit_finalO18rate01(:,3)),L,'false');
    V_01_pick=Orbit_finalO18rate01(ID_picked,3);
    V_02_pick=Orbit_finalO18rate02(ID_picked,3);
    Rand_V_01_pick=V_01_pick(randsample(L,L,'false'));
    Rand_V_02_pick=V_02_pick(randsample(L,L,'false'));
    Rand_Orbit_o18rate01=Orbit_finalO18rate01(:,3);
    Rand_Orbit_o18rate02=Orbit_finalO18rate02(:,3);
    Rand_Orbit_o18rate01(ID_picked)=Rand_V_01_pick;
    Rand_Orbit_o18rate02(ID_picked)=Rand_V_02_pick;
    STD_rand(num)=std(log2(Rand_Orbit_o18rate01)-log2(Rand_Orbit_o18rate02));
    MEAN_rand(num)=mean(log2(Rand_Orbit_o18rate01)-log2(Rand_Orbit_o18rate02));
    num=num+1;
end
figure
plot(1:20:200,STD_rand);grid on;
figure
plot(1:20:200,MEAN_rand);grid on;
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%plot 3D spectrum
for i=15;
    ID_into=Com_ID_0102(Training_id(i));
    Sp_mass=TOF_mass_com(Training_id(i));
    Sp_cs=TOF_charge_com(Training_id(i));
    data01_intervalid=Training_Row_Column_Num(i,1);
    data02_intervalid=Training_Row_Column_Num(i,2);
    TOF_tr_scan_start01=IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,1);
    TOF_tr_scan_end01=IntervalList01{ID_into}.intervallist_after7_combine(data01_intervalid,2);
    TOF_tr_scan_start02=IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,1);
    TOF_tr_scan_end02=IntervalList02{ID_into}.intervallist_after7_combine(data02_intervalid,2);
    mzvalue=(Sp_mass+Sp_cs*1.0073)/Sp_cs;
    MZ_interval=[mzvalue-mzvalue*80/10^6 (mzvalue+5.01184/Sp_cs)+(mzvalue+5.01184/Sp_cs)*80/10^6];
    Interval_start=1500;
    Interval_end=2000;
%     Interval_start=TOF_tr_scan_start01;
%     Interval_end=TOF_tr_scan_end01;
    [Int_Matrix_mz_time,mz_unique,Pep_mz,Interval_XIC]=Generate_Pep_Mz_Time_matrix(TOF_peakl01,MZ_interval,Interval_start,Interval_end);
    MZ_AXIS=mz_unique*ones(1,size(Int_Matrix_mz_time,2));
    T_unique=TOF_retentiont01l1(Interval_start:Interval_end);
    T_AXIS=ones(size(Int_Matrix_mz_time,1),1)*T_unique;
    figure
    surf(MZ_AXIS,T_AXIS,Int_Matrix_mz_time)
    figure
    plot(mz_unique,sum(Int_Matrix_mz_time,2))
    grid on
    
%     Interval_start=TOF_tr_scan_start02;
%     Interval_end=TOF_tr_scan_end02;
    [Int_Matrix_mz_time,mz_unique,Pep_mz,Interval_XIC]=Generate_Pep_Mz_Time_matrix(TOF_peakl02,MZ_interval,Interval_start,Interval_end);
    MZ_AXIS=mz_unique*ones(1,size(Int_Matrix_mz_time,2));
    T_unique=TOF_retentiont02l1(Interval_start:Interval_end);
    T_AXIS=ones(size(Int_Matrix_mz_time,1),1)*T_unique;
   
    figure
    surfc(MZ_AXIS,T_AXIS,Int_Matrix_mz_time)
    figure
    plot(mz_unique,sum(Int_Matrix_mz_time,2))
    grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Times_noise_std=3;%6;
instrument_h_int=10^7;
de_range=10^5;%10^4;
colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
for i=151:155%%length(posi01v2)
    data01_intervalid=Training_Row_Column_Num(i,1);
    data02_intervalid=Training_Row_Column_Num(i,2);
    ID_into=Com_ID_0102(Training_id(i));
    Labeling_eff_diff=abs(IntervalList01{ID_into}.Labelling_efficiency_after7_combine-Orbit_Labelling_efficiency(ID_into));
   
%     for j=1:6
%         th(j)=TOF_getThreshold_range(TOF_XICs01{j}(:,ID_into),instrument_h_int, de_range,Times_noise_std);
%         intervalRawData(:,j)=TOF_XICs01{j}(IntervalList01{ID_into}.intervallist_after7_combine(2,1):IntervalList01{ID_into}.intervallist_after7_combine(2,2),ID_into);
%     end
%     Judge=long_criteria(iso(ID_into,:),intervalRawData,th);
    
    figure
    subplot(2,2,1)
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
    for kkk=1:length(IntervalList01{ID_into}.intervallist_after7_combine(:,1))
        interval_mid_time=(TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(kkk,1))+TOF_retentiont01l1(IntervalList01{ID_into}.intervallist_after7_combine(kkk,2)))/2;
        text(interval_mid_time,2*Height,num2str(Labeling_eff_diff(kkk)));
    end
    grid on
    subplot(2,2,2)
    for j=1:6
        plot(retentiont01l1,Orbit_TotalXICs01{j}(:,ID_into),colorarray(j))
        hold on
    end
    stem(retentiont01l1(Orbit_AMTms2time(ID_into,2)),10^7,'ro')
    stem(retentiont01l1(Orbit_AMTms2time(ID_into,3)),10^7,'ko')
    grid on    
    
    clear Labeling_eff_diff
    Labeling_eff_diff=abs(IntervalList02{ID_into}.Labelling_efficiency_after7_combine-Orbit_Labelling_efficiency(ID_into));
    subplot(2,2,3)
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
    for kkk=1:length(IntervalList02{ID_into}.intervallist_after7_combine(:,1))
        interval_mid_time=(TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(kkk,1))+TOF_retentiont02l1(IntervalList02{ID_into}.intervallist_after7_combine(kkk,2)))/2;
        text(interval_mid_time,2*Height,num2str(Labeling_eff_diff(kkk)));
    end
    grid on
    
    clear Labeling_eff_diff
    
%     MZ_interval(1)=TOF_totalmzList(ID_into,1)-0.2;
%     MZ_interval(2)=TOF_totalmzList(ID_into,6)+0.2;
%     Interval_start=IntervalList01{ID_into}.intervallist_after7_combine(1,1);
%     Interval_end=IntervalList01{ID_into}.intervallist_after7_combine(1,2);
%     [Int_Matrix_mz_time,mz_unique,Pep_mz,Interval_XIC]=Generate_Pep_Mz_Time_matrix(TOF_peakl01,MZ_interval,Interval_start,Interval_end);
%     x=mz_unique*ones(1,size(Int_Matrix_mz_time,2));
%     y=ones(size(Int_Matrix_mz_time,1),1)*[Interval_start:Interval_end];
%     figure
%     surf(x,y,Int_Matrix_mz_time)    
%     figure
%     plot(mz_unique,sum(Int_Matrix_mz_time,2))
%     hold on
%     grid on
%     stem(TOF_totalmzList(ID_into,:),30000*ones(1,6),'ro')
    
end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% based on training data generate R T LE model for
%%%%%%%%%%%%%%%%%% TOF2Orbit
for i=1:length(Training_id)
    Training_TOFInterval_id01(i)=Training_Row_Column_Num(i,1);
    k=Training_id(i);
    Training_TOFscore01{i}.LOGKL_TOFOrbit=TOFscorev1{k}.LOGKL_TOF01Orbit01;
    Training_TOFscore01{i}.R_TOFOrbit=TOFscorev1{k}.R_TOF01Orbit01;
    Training_TOFscore01{i}.T_TOFOrbit=TOFscorev1{k}.T_TOF01Orbit01;
    Training_TOFscore01{i}.Normal_T_TOFOrbit=TOFscorev1{k}.Normal_T_TOF01Orbit01;
    Training_TOFscore01{i}.T_TOF=TOFscorev1{k}.T_TOF01;
    Training_TOFscore01{i}.Normal_T_TOF=TOFscorev1{k}.Normal_T_TOF01;
    Training_TOFscore01{i}.T_Orbit=TOFscorev1{k}.T_Orbit01;
    Training_TOFscore01{i}.Normal_T_Orbit=TOFscorev1{k}.Normal_T_Orbit01;
    Training_TOFscore01{i}.LE_TOF=TOFscorev1{k}.LE_TOF01;
    Training_TOFscore01{i}.LE_TOFOrbit=TOFscorev1{k}.LE_TOF01Orbit01;    
    Training_Intervallist01(i)=IntervalList_com01(k);
    Training_Orbitms2time01(i)=TOF_ms2information_com(k);
end

[T_null_TOF01_Orbit01,R_null_TOF01_Orbit01,LE_null_TOF01_Orbit01...
    PP_null_TOF2Orbit,PP_null_Normal_TOF2Orbit,...
    Mu_null_TOF2Orbit,Sigma_null_TOF2Orbit,...
    PHAT_null,T_null_95per_Bond,R_null_95per_Bond]=Get_TOF_Orbit_RTLE_null_model(Training_Orbitms2time01,...
            Training_Intervallist01,retentiont01l1,Training_TOFscore01,Std_AMT_timeshift,Training_TOFInterval_id01);



%%%%%%%%%%%%%%%%%% Get R T LE null model and pep id based on TOF01 and
%%%%%%%%%%%%%%%%%% Orbit01
%%%%%%%%%%%%%%%%%%
for i=1:length(TOFscorev1)
    TOFscore01{i}.LOGKL_TOFOrbit=TOFscorev1{i}.LOGKL_TOF01Orbit01;
    TOFscore01{i}.R_TOFOrbit=TOFscorev1{i}.R_TOF01Orbit01;
    TOFscore01{i}.T_TOFOrbit=TOFscorev1{i}.T_TOF01Orbit01;
    TOFscore01{i}.Normal_T_TOFOrbit=TOFscorev1{i}.Normal_T_TOF01Orbit01;
    TOFscore01{i}.T_TOF=TOFscorev1{i}.T_TOF01;
    TOFscore01{i}.Normal_T_TOF=TOFscorev1{i}.Normal_T_TOF01;
    TOFscore01{i}.T_Orbit=TOFscorev1{i}.T_Orbit01;
    TOFscore01{i}.Normal_T_Orbit=TOFscorev1{i}.Normal_T_Orbit01;
    TOFscore01{i}.LE_TOF=TOFscorev1{i}.LE_TOF01;
    TOFscore01{i}.LE_TOFOrbit=TOFscorev1{i}.LE_TOF01Orbit01;   
end

[T_null_TOF01_Orbit01,R_null_TOF01_Orbit01,PP_TOF2Orbit01,PP_Normal_TOF2Orbit01,Mu_TOF2Orbit01,Sigma_TOF2Orbit01,R_PHAT_TOF2Orbit01,T_null_95per_Bond_TOF2Orbit01,R_null_95per_Bond_TOF2Orbit01,Null_model_id_TOF2Orbit01]=Get_TOF_Orbit_RT_null_model(Pep_com,TOF_ms2information_com,IntervalList_com01,retentiont01l1,TOFscore01,Std_AMT_timeshift);

% [R_Training_ID01,T_TOF01_Orbit01,T_Training_ID01,R_TOF01_Orbit01,PP_TOF2Orbit01,PP_Normal_TOF2Orbit01,Mu_TOF2Orbit01,Sigma_TOF2Orbit01]=Get_TOF_Orbit_RT_model(Pep_com,TOF_ms2information_com,IntervalList_com01,retentiont01l1,TOFscore01,Std_AMT_timeshift);

b1=-1:0.01:1;
Y1_sig=normpdf(b1,Mu_TOF2Orbit01,Sigma_TOF2Orbit01);
figure
plot(b1,Y1_sig);
b1=0:0.05:1;
Y1_sig=gampdf(1-b1,R_PHAT_TOF2Orbit01(1),R_PHAT_TOF2Orbit01(2));
figure
plot(b1,Y1_sig);
Training_id_R_T_01=Null_model_id_TOF2Orbit01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Get R_Training_data based on TOF02 and Orbit02
%%%%%%%%%%%%%%%%%%
for i=1:length(TOFscorev1)
    TOFscore02{i}.LOGKL_TOFOrbit=TOFscorev1{i}.LOGKL_TOF02Orbit02;
    TOFscore02{i}.R_TOFOrbit=TOFscorev1{i}.R_TOF02Orbit02;
    TOFscore02{i}.T_TOFOrbit=TOFscorev1{i}.T_TOF02Orbit02;
    TOFscore02{i}.Normal_T_TOFOrbit=TOFscorev1{i}.Normal_T_TOF02Orbit02;
    TOFscore02{i}.T_TOF=TOFscorev1{i}.T_TOF02;
    TOFscore02{i}.Normal_T_TOF=TOFscorev1{i}.Normal_T_TOF02;
    TOFscore02{i}.T_Orbit=TOFscorev1{i}.T_Orbit02;
    TOFscore02{i}.Normal_T_Orbit=TOFscorev1{i}.Normal_T_Orbit02;
end
[T_null_TOF02_Orbit02,R_null_TOF02_Orbit02,PP_TOF2Orbit02,PP_Normal_TOF2Orbit02,Mu_TOF2Orbit02,Sigma_TOF2Orbit02,R_PHAT_TOF2Orbit02,T_null_95per_Bond_TOF2Orbit02,R_null_95per_Bond_TOF2Orbit02,Null_model_id_TOF2Orbit02]=Get_TOF_Orbit_RT_null_model(Pep_com,TOF_ms2information_com,IntervalList_com02,retentiont02l1,TOFscore02,Std_AMT_timeshift);

% [R_Training_ID01,T_TOF01_Orbit01,T_Training_ID01,R_TOF01_Orbit01,PP_TOF2Orbit01,PP_Normal_TOF2Orbit01,Mu_TOF2Orbit01,Sigma_TOF2Orbit01]=Get_TOF_Orbit_RT_model(Pep_com,TOF_ms2information_com,IntervalList_com01,retentiont01l1,TOFscore01,Std_AMT_timeshift);

b1=-1:0.01:1;
Y1_sig=normpdf(b1,Mu_TOF2Orbit02,Sigma_TOF2Orbit02);
figure
plot(b1,Y1_sig);
b1=0:0.05:1;
Y1_sig=gampdf(1-b1,R_PHAT_TOF2Orbit02(1),R_PHAT_TOF2Orbit02(2));
figure
plot(b1,Y1_sig);
Training_id_R_T_02=Null_model_id_TOF2Orbit02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length(Training_id_R_T_01)
length(Training_id_R_T_02)
% length(intersect(Training_id_R_T_01,Training_id_R_T_02))
Training_id_R_T=intersect(Training_id_R_T_01,Training_id_R_T_02);
length(Training_id_R_T)


% %%%%%%%%%%%%%%%%%% generate AR AT AKL model
% %%%%%%%%%%%%%%%%%% select training data only based on QTOF data
% Timepair_corr_basedonAR=[];
% Timepair_corr_basedonAR_normal=[];
% Timepair_corr_basedonAKL=[];
% Timepair_corr_basedonAKL_normal=[];
% Timepair_noncorr_basedonAR=[];
% Timepair_noncorr_basedonAR_normal=[];
% Timepair_noncorr_basedonAKL=[];
% Timepair_noncorr_basedonAKL_normal=[];
% % AT_corr=[];AT_non_corr=[];
% AR_corr=[];AR_non_corr=[];
% AKL0102_corr=[];AKL0102_non_corr=[];
% Training_id=[];
% 
% for i=1:length(TOFscore)
%     Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
%     Interval_matrix02=IntervalList_com02{i}.intervallist_after7_combine;
%     
%     [R_AT,C_AT,V_AT]=find(abs(TOFscore{i}.Normal_AT)<=Std_AMT_timeshift);
%     if ~isempty(R_AT)
%         V_LOGAKL=diag(TOFscore{i}.LOGAKL(R_AT,C_AT));
%         V_AR=diag(TOFscore{i}.AR(R_AT,C_AT));
%         AKL_good_id=find(V_LOGAKL<=-5);
%         AR_good_id=find(V_AR>=0.95);
%         Intersection_AKL_AR=intersect(AKL_good_id,AR_good_id);
%         
%         if ~isempty(Intersection_AKL_AR)%~isempty(AKL_good_id) && ~isempty(AR_good_id)
% %             R_AT_AKL_AR=R_AT(Intersection_AKL_AR);
% %             C_AT_AKL_AR=C_AT(Intersection_AKL_AR);
% %             for k=1:length(R_AT_AKL_AR)
% %                 Sc_Start01=Interval_matrix01(R_AT_AKL_AR(k),1);
% %                 Sc_end01=Interval_matrix01(R_AT_AKL_AR(k),2);
% %                 Sc_Start02=Interval_matrix02(C_AT_AKL_AR(k),1);
% %                 Sc_end02=Interval_matrix02(C_AT_AKL_AR(k),2);
% %             end
%             
%                 Training_id=[Training_id;i];
%                     %%%%%%% AR corr and non-corr
%                     [V_LOGAKL_sort,ID_AKL_sort]=sort(V_LOGAKL);
%                     UP_b=ceil(length(V_LOGAKL)/3);
%                     [V_AR_max,P_AR_max]=max(V_AR(ID_AKL_sort(1:UP_b)));
%                     R_corr=R_AT(ID_AKL_sort(P_AR_max));
%                     C_corr=C_AT(ID_AKL_sort(P_AR_max));
%                     Interval_time01=sum(TOF_retentiont01l1(Interval_matrix01),2)/2;
%                     Interval_time02=sum(TOF_retentiont02l1(Interval_matrix02),2)/2;
%                     Timepoint01=Interval_time01(R_corr);%sum(TOF_retentiont01l1(Interval_matrix01(R_corr,:)))/2;
%                     Timepoint02=Interval_time02(C_corr);%sum(TOF_retentiont02l1(Interval_matrix02(C_corr,:)))/2;
% %                     Timepoint01_noncorr=Interval_time01;Timepoint01_noncorr(R_corr)=[];
%                     Timepoint02_noncorr=Interval_time02;Timepoint02_noncorr(C_corr)=[];
%                     Timepoint01_normal=Timepoint01/max(TOF_retentiont01l1);
%                     Timepoint02_normal=Timepoint02/max(TOF_retentiont02l1);
%                     Timepair_corr_basedonAR=[Timepair_corr_basedonAR;Timepoint01 Timepoint02];
%                     Timepair_corr_basedonAR_normal=[Timepair_corr_basedonAR_normal; Timepoint01_normal Timepoint02_normal];
%                     Timepair_noncorr_basedonAR=[Timepair_noncorr_basedonAR;Timepoint01.*ones(length(Timepoint02_noncorr),1) Timepoint02_noncorr];
%                     Timepair_noncorr_basedonAR_normal=[Timepair_noncorr_basedonAR_normal;Timepoint01.*ones(length(Timepoint02_noncorr),1)./max(TOF_retentiont01l1) Timepoint02_noncorr./max(TOF_retentiont02l1)];                    
%                     Vector_AR_noncorr_Row=TOFscore{i}.AR(R_corr,:);
%                     Vector_AR_noncorr_Col=TOFscore{i}.AR(:,C_corr);
%                     Vector_AR_noncorr_Row(C_corr)=[];
%                     Vector_AR_noncorr_Col(R_corr)=[];
%                     AR_corr=[AR_corr;V_AR_max];
%                     AR_non_corr=[AR_non_corr;Vector_AR_noncorr_Row';Vector_AR_noncorr_Col];
%                     %%%%%%%
% 
%                     %%%%%%% AKL corr and non-corr
%                     [V_AR_sort,ID_AR_sort]=sort(V_AR,'descend');
%                     UP_b=ceil(length(V_AR)/3);
%                     [V_AKL_max,P_AKL_max]=min(V_LOGAKL(ID_AR_sort(1:UP_b)));
%                     R_corr=R_AT(ID_AR_sort(P_AKL_max));
%                     C_corr=C_AT(ID_AR_sort(P_AKL_max));
%                     Interval_time01=sum(TOF_retentiont01l1(Interval_matrix01),2)/2;
%                     Interval_time02=sum(TOF_retentiont02l1(Interval_matrix02),2)/2;
%                     Timepoint01=Interval_time01(R_corr);%sum(TOF_retentiont01l1(Interval_matrix01(R_corr,:)))/2;
%                     Timepoint02=Interval_time02(C_corr);%sum(TOF_retentiont02l1(Interval_matrix02(C_corr,:)))/2;
%                     Timepoint02_noncorr=Interval_time02;Timepoint02_noncorr(C_corr)=[];
%                     Timepoint01_normal=Timepoint01/max(TOF_retentiont01l1);
%                     Timepoint02_normal=Timepoint02/max(TOF_retentiont02l1);
%                     Timepair_corr_basedonAKL=[Timepair_corr_basedonAKL;Timepoint01 Timepoint02];
%                     Timepair_corr_basedonAKL_normal=[Timepair_corr_basedonAKL_normal; Timepoint01_normal Timepoint02_normal];
%                     Timepair_noncorr_basedonAKL=[Timepair_noncorr_basedonAKL;Timepoint01.*ones(length(Timepoint02_noncorr),1) Timepoint02_noncorr];
%                     Timepair_noncorr_basedonAKL_normal=[Timepair_noncorr_basedonAKL_normal;Timepoint01.*ones(length(Timepoint02_noncorr),1)./max(TOF_retentiont01l1) Timepoint02_noncorr./max(TOF_retentiont02l1)];                    
%                     Vector_AKL_noncorr_Row=TOFscore{i}.LOGAKL(R_corr,:);
%                     Vector_AKL_noncorr_Col=TOFscore{i}.LOGAKL(:,C_corr);
%                     Vector_AKL_noncorr_Row(C_corr)=[];
%                     Vector_AKL_noncorr_Col(R_corr)=[];
%                     AKL0102_corr=[AKL0102_corr;V_AKL_max];
%                     AKL0102_non_corr=[AKL0102_non_corr;Vector_AKL_noncorr_Row';Vector_AKL_noncorr_Col];        
%                     %%%%%%%
%         end
%     end
% end
% 
% % save Training_id_basedon_QTOF Training_id
% %%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%% select training data based on Orbit ms2 time
% %%%%%%%%%%%%%%%%%%%%%%% information
% % Timepair_corr_basedonAR=[];
% % Timepair_corr_basedonAR_normal=[];
% % Timepair_corr_basedonAKL=[];
% % Timepair_corr_basedonAKL_normal=[];
% % Timepair_noncorr_basedonAR=[];
% % Timepair_noncorr_basedonAR_normal=[];
% % Timepair_noncorr_basedonAKL=[];
% % Timepair_noncorr_basedonAKL_normal=[];
% Timepair_corr=[];
% Timepair_corr_normal=[];
% Timepair_noncorr=[];
% Timepair_noncorr_normal=[];
% AT_corr=[];AT_non_corr=[];
% AR_corr=[];AR_non_corr=[];
% AKL0102_corr=[];AKL0102_non_corr=[];
% Training_id=[];
% for i=1:length(TOFscore)
%     Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
%     Interval_matrix02=IntervalList_com02{i}.intervallist_after7_combine;
%          
%     Normalized_ms2time01_TOF=TOF_ms2information_com(i)/max(retentiont01l1);
%     Normalized_ms2time02_TOF=TOF_ms2information_com(i)/max(retentiont02l1);
%     
%     R_AT=find(TOFscore{i}.Normal_T01(:,1)<=Normalized_ms2time01_TOF+1*Std_AMT_timeshift & TOFscore{i}.Normal_T01(:,1)>=Normalized_ms2time01_TOF-1*Std_AMT_timeshift);
%     C_AT=find(TOFscore{i}.Normal_T02(1,:)<=Normalized_ms2time02_TOF+1*Std_AMT_timeshift & TOFscore{i}.Normal_T02(1,:)>=Normalized_ms2time02_TOF-1*Std_AMT_timeshift);
%     %%%%% 1. Find any intervals in TOF that is in Orbit ms2 interval
%     %%%%% 2. LOGAKL <= -2.5
%     %%%%% 3. AR >= 0.85
% %     if ~isempty(R_AT) && ~isempty(C_AT)
%     %%%%% 1. Orbit ms2 interval only contains one interval in TOF
%     %%%%% 2. LOGAKL <= -2.5
%     %%%%% 3. AR >= 0.85
%     if length(R_AT)==1 && length(C_AT)==1
%         
%         V_LOGAKL=TOFscore{i}.LOGAKL(R_AT,C_AT);
%         V_AR=TOFscore{i}.AR(R_AT,C_AT);
%         
%         [AKL_R_good_id,AKL_C_good_id,AKL_V_good_id]=find(V_LOGAKL<=-2.5);
%         [AR_R_good_id,AR_C_good_id,AR_V_good_id]=find(V_AR>=0.85);
%         
%         AKL_good_id_vector=[AKL_R_good_id,AKL_C_good_id];
%         AR_good_id_vector=[AR_R_good_id,AR_C_good_id];
%         
%         Good_R_id=[];
%         Good_C_id=[];
%         if ~isempty(AKL_good_id_vector) && ~isempty(AR_good_id_vector)
%             for j=1:size(AKL_good_id_vector,1)
%                 index_diff01=AR_good_id_vector(:,1)-AKL_good_id_vector(j,1);
%                 index_diff02=AR_good_id_vector(:,2)-AKL_good_id_vector(j,2);
%                 Same_index=find(abs(index_diff01)+abs(index_diff02)==0);
%                 if ~isempty(Same_index)
%                     Good_R_id=[Good_R_id;AKL_good_id_vector(j,1)];
%                     Good_C_id=[Good_C_id;AKL_good_id_vector(j,2)];
%                 end
%             end
%             if ~isempty(Good_R_id)
%                 R_corr=R_AT(Good_R_id);
%                 C_corr=C_AT(Good_C_id);
%                 if length(R_corr)==1 && length(C_corr)==1
% 
%                     Training_id=[Training_id; i];
% 
%                     Vector_T01_noncorr_Row=TOFscore{i}.T01(:,1);
%                     Vector_T02_noncorr_Col=TOFscore{i}.T02(1,:);
%                     Vector_T01_noncorr_Row(R_corr)=[];
%                     Vector_T02_noncorr_Col(C_corr)=[];
%                     Vector_T01_normal_noncorr_Row=TOFscore{i}.Normal_T01(:,1);
%                     Vector_T02_normal_noncorr_Col=TOFscore{i}.Normal_T02(1,:);
%                     Vector_T01_normal_noncorr_Row(R_corr)=[];
%                     Vector_T02_normal_noncorr_Col(C_corr)=[];
%                     Timepair_corr=[Timepair_corr;TOFscore{i}.T01(R_corr,C_corr),TOFscore{i}.T02(R_corr,C_corr)];
%                     Timepair_corr_normal=[Timepair_corr_normal;TOFscore{i}.Normal_T01(R_corr,C_corr),TOFscore{i}.Normal_T02(R_corr,C_corr)];
%                     Timepair_noncorr=[Timepair_noncorr;TOFscore{i}.T01(R_corr,C_corr)*ones(length(Vector_T02_noncorr_Col),1),Vector_T02_noncorr_Col'];
%                     Timepair_noncorr_normal=[Timepair_noncorr_normal;TOFscore{i}.Normal_T01(R_corr,C_corr)*ones(length(Vector_T02_normal_noncorr_Col),1),Vector_T02_normal_noncorr_Col'];
%                     
%                     
%                     V_AT_corr=TOFscore{i}.Normal_AT(R_corr,C_corr);
%                     Vector_AT_noncorr_Row=TOFscore{i}.Normal_AT(R_corr,:);
%                     Vector_AT_noncorr_Col=TOFscore{i}.Normal_AT(:,C_corr);
%                     Vector_AT_noncorr_Row(C_corr)=[];
%                     Vector_AT_noncorr_Col(R_corr)=[];
%                     AT_corr=[AT_corr;V_AT_corr];
%                     AT_non_corr=[AT_non_corr;Vector_AT_noncorr_Row';Vector_AT_noncorr_Col];
% 
%                     V_AR_corr=TOFscore{i}.AR(R_corr,C_corr);
%                     Vector_AR_noncorr_Row=TOFscore{i}.AR(R_corr,:);
%                     Vector_AR_noncorr_Col=TOFscore{i}.AR(:,C_corr);
%                     Vector_AR_noncorr_Row(C_corr)=[];
%                     Vector_AR_noncorr_Col(R_corr)=[];
%                     AR_corr=[AR_corr;V_AR_corr];
%                     AR_non_corr=[AR_non_corr;Vector_AR_noncorr_Row';Vector_AR_noncorr_Col];
% 
%                     V_AKL_corr=TOFscore{i}.LOGAKL(R_corr,C_corr);
%                     Vector_AKL_noncorr_Row=TOFscore{i}.LOGAKL(R_corr,:);
%                     Vector_AKL_noncorr_Col=TOFscore{i}.LOGAKL(:,C_corr);
%                     Vector_AKL_noncorr_Row(C_corr)=[];
%                     Vector_AKL_noncorr_Col(R_corr)=[];
%                     AKL0102_corr=[AKL0102_corr;V_AKL_corr];
%                     AKL0102_non_corr=[AKL0102_non_corr;Vector_AKL_noncorr_Row';Vector_AKL_noncorr_Col];        
%                 end
%             end
%         end
%     end    
% end
% length(Training_id)
% save Training_id_basedon_Orbit  Training_id
% %%%%%%%%%%%%%%%%%%%%%%%
figure
plot(Timepair_corr(:,1),Timepair_corr(:,2),'r*')

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

[AKL0102_MU_sig,AKL0102_SIGMA_sig]=normfit(AKL0102_corr);
[AKL0102_MU_notsig,AKL0102_SIGMA_notsig]=normfit(AKL0102_non_corr);
xx=-20:0.1:10;
figure
subplot(2,2,1)
hist(AKL0102_corr,xx);grid on;
subplot(2,2,3)
hist(AKL0102_non_corr,xx);grid on;
subplot(1,2,2)
plot(xx,normpdf(xx,AKL0102_MU_sig,AKL0102_SIGMA_sig))
hold on
plot(xx,normpdf(xx,AKL0102_MU_notsig,AKL0102_SIGMA_notsig),'r');grid on;


PP=polyfit(Timepair_corr_normal(:,1),Timepair_corr_normal(:,2),4);
figure
plot(Timepair_corr_normal(:,1),Timepair_corr_normal(:,2),'r.')
XXX=0:0.01:1;
YYY=polyval(PP,XXX);
hold on
plot(XXX,YYY)

% Timepair_corr_basedonAR;
AT_corr=polyval(PP,Timepair_corr_normal(:,1))-Timepair_corr_normal(:,2);
% Timepair_noncorr_basedonAR;
AT_non_corr=polyval(PP,Timepair_noncorr_normal(:,1))-Timepair_noncorr_normal(:,2);
[timemodel0102_MU_sig,timemodel0102_SIGMA_sig] = normfit(AT_corr);
[timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig] = normfit(AT_non_corr);
t=-1:0.001:1;
P_timenorm_sigv2=normpdf(t,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
P_timenorm_notsigv2=normpdf(t,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
figure
plot(t,P_timenorm_sigv2,'r');hold on
plot(t,P_timenorm_notsigv2)

%%%%%%%%%%%%%%%%%% Generate Testing data

Testing_id=1:length(TOFscorev1);
Testing_id(Training_id)=[];
Testingdata=TOFscorev1(Testing_id);

iso_training=iso_com(Training_id);
Training_Orbit_AMTms2time_com=Orbit_AMTms2time_com(Training_id);
for i=1:6
    Orbit_XICs_training01{i}=Orbit_XICs_com01{i}(:,Training_id);
    TOF_XICs_training01{i}=TOF_XICs_com01{i}(:,Training_id);
end
IntervalList_training01=IntervalList_com01(Training_id);
for i=1:6
    Orbit_XICs_training02{i}=Orbit_XICs_com02{i}(:,Training_id);
    TOF_XICs_training02{i}=TOF_XICs_com02{i}(:,Training_id);
end
IntervalList_training02=IntervalList_com02(Training_id);

iso_testing=iso_com(Testing_id,:);
Testing_Orbit_AMTms2time_com=Orbit_AMTms2time_com(Testing_id,:);
for i=1:6
    Orbit_XICs_testing01{i}=Orbit_XICs_com01{i}(:,Testing_id);
    TOF_XICs_testing01{i}=TOF_XICs_com01{i}(:,Testing_id);
end
IntervalList_testing01=IntervalList_com01(Testing_id);
for i=1:6
    Orbit_XICs_testing02{i}=Orbit_XICs_com02{i}(:,Testing_id);
    TOF_XICs_testing02{i}=TOF_XICs_com02{i}(:,Testing_id);
end
IntervalList_testing02=IntervalList_com02(Testing_id);
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
num=0;
num_notdetect=0;
Orbit_goodalign_on_testing=[];
TOF_goodalign_corrBnon_on_testing=[];
TOF_Orbit_goodalign0102_timescale_on_testing=[];
TOF_Orbit_goodalign01_on_testing=[];
TOF_Orbit_goodalign02_on_testing=[];
TOF_Orbit_goodalign0102_on_testing=[];

for i=1:length(Testingdata)
    
    Orbit01_ms2_interval_start=Testing_Orbit_AMTms2time_com(i,2);
    Orbit01_ms2_interval_end=Testing_Orbit_AMTms2time_com(i,3);
    Orbit02_ms2_interval_start=Testing_Orbit_AMTms2time_com(i,5);
    Orbit02_ms2_interval_end=Testing_Orbit_AMTms2time_com(i,6);

    if Orbit01_ms2_interval_start~=0 && Orbit02_ms2_interval_start~=0
        Orbit01_ms2_time=(retentiont01l1(Orbit01_ms2_interval_start)+retentiont01l1(Orbit01_ms2_interval_end))/2;
        Orbit02_ms2_time=(retentiont02l1(Orbit02_ms2_interval_start)+retentiont02l1(Orbit02_ms2_interval_end))/2;

        Orbit_goodalign_on_testing=[Orbit_goodalign_on_testing;i];
        Pepinformdetect(i).Orbit_align_JUD=1;        

        AR_corr_Prob=betapdf(Testingdata{i}.AR,AR0102_PARMHAT1(1),AR0102_PARMHAT1(2));
        AR_non_corr_Prob=betapdf(Testingdata{i}.AR,AR0102_PARMHAT2(1),AR0102_PARMHAT2(2));
        AKL_corr_Prob=normpdf(Testingdata{i}.LOGAKL,AKL0102_MU_sig,AKL0102_SIGMA_sig);
        AKL_non_corr_Prob=normpdf(Testingdata{i}.LOGAKL,AKL0102_MU_notsig,AKL0102_SIGMA_notsig);
        AT_corr_Prob=normpdf(polyval(PP,Testingdata{i}.Normal_T01)-Testingdata{i}.Normal_T02,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
        AT_non_corr_Prob=normpdf(polyval(PP,Testingdata{i}.Normal_T01)-Testingdata{i}.Normal_T02,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
        Prob_AR_matrix=log(AR_corr_Prob./AR_non_corr_Prob);
        Prob_AKL_matrix=log(AKL_corr_Prob./AKL_non_corr_Prob);
        Prob_AT_matrix=log(AT_corr_Prob./AT_non_corr_Prob);
        Total_Prob_matrix=Prob_AR_matrix+Prob_AKL_matrix+Prob_AT_matrix;
        [R_detect,C_detect,V_detect]=find(Total_Prob_matrix>=0);
        [V_max_detect,I_max_detect]=max(Total_Prob_matrix);
        [V_max_detectv1,I_max_detectv1]=max(V_max_detect);    
        R_max_detect=I_max_detect(I_max_detectv1);
        C_max_detect=I_max_detectv1;
        if ~isempty(R_detect)
%             num=num+1;
            TOF_goodalign_corrBnon_on_testing=[TOF_goodalign_corrBnon_on_testing;i];
            Pepinformdetect(i).ID_detect=i;
            Pepinformdetect(i).R_detect=R_detect;
            Pepinformdetect(i).C_detect=C_detect;
            Pepinformdetect(i).R_max_detect=R_max_detect;
            Pepinformdetect(i).C_max_detect=C_max_detect;
            
            for mzid=1:6
                Orbit_XICs_matrix01(:,mzid)=Orbit_XICs_testing01{mzid}(Orbit01_ms2_interval_start:Orbit01_ms2_interval_end,i);
                Orbit_XICs_matrix02(:,mzid)=Orbit_XICs_testing02{mzid}(Orbit02_ms2_interval_start:Orbit02_ms2_interval_end,i);
            end
            ID1=0; ID2=0; ID3=0; ID4=0;
            for k=1:length(R_detect)
                TOF01_start=IntervalList_testing01{i}.intervallist_after7_combine(R_detect(k),1);
                TOF01_end=IntervalList_testing01{i}.intervallist_after7_combine(R_detect(k),2);
                TOF02_start=IntervalList_testing02{i}.intervallist_after7_combine(C_detect(k),1);
                TOF02_end=IntervalList_testing02{i}.intervallist_after7_combine(C_detect(k),2);
                TOF01_start_time=TOF_retentiont01l1(TOF01_start);
                TOF01_end_time=TOF_retentiont01l1(TOF01_end);
                TOF02_start_time=TOF_retentiont02l1(TOF02_start);
                TOF02_end_time=TOF_retentiont02l1(TOF02_end);
                if polyval(PP_Normal_TOF2Orbit01,(sum(TOF01_start_time+TOF01_end_time)/2)/max(TOF_retentiont01l1))<=Orbit01_ms2_time/max(retentiont01l1)+1*Std_AMT_timeshift && polyval(PP_Normal_TOF2Orbit01,(sum(TOF01_start_time+TOF01_end_time)/2)/max(TOF_retentiont01l1))>=Orbit01_ms2_time/max(retentiont01l1)-1*Std_AMT_timeshift
                    if polyval(PP_Normal_TOF2Orbit02,(sum(TOF02_start_time+TOF02_end_time)/2)/max(TOF_retentiont02l1))<=Orbit02_ms2_time/max(retentiont02l1)+1*Std_AMT_timeshift && polyval(PP_Normal_TOF2Orbit02,(sum(TOF02_start_time+TOF02_end_time)/2)/max(TOF_retentiont02l1))>=Orbit02_ms2_time/max(retentiont02l1)-1*Std_AMT_timeshift
                        ID4=ID4+1;
                        TOF_Orbit_goodalign0102_timescale_on_testing=[TOF_Orbit_goodalign0102_timescale_on_testing;i];
                        Pepinformdetect(i).R_TOF_Orbit_timescale_detect(ID4)=R_detect(k);
                        Pepinformdetect(i).C_TOF_Orbit_timescale_detect(ID4)=C_detect(k);
                            for mzid=1:6
                                TOF_XICs_matrix01(:,mzid)=TOF_XICs_testing01{mzid}(TOF01_start:TOF01_end,i);
                                TOF_XICs_matrix02(:,mzid)=TOF_XICs_testing02{mzid}(TOF02_start:TOF02_end,i);
                                if sum(TOF_XICs_matrix01(:,mzid))==0 || sum(Orbit_XICs_matrix01(:,mzid))==0
                                    AR_TOF_Orbit_01(mzid)=0;
                                else
                                    [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( TOF_XICs_matrix01(:,mzid)', Orbit_XICs_matrix01(:,mzid)');
                                    [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                                    AR_TOF_Orbit_01(mzid)=STATS(1);                        
                                end
                                if sum(TOF_XICs_matrix02(:,mzid))==0 || sum(Orbit_XICs_matrix02(:,mzid))==0
                                    AR_TOF_Orbit_02(mzid)=0;
                                else
                                    [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( TOF_XICs_matrix02(:,mzid)', Orbit_XICs_matrix02(:,mzid)');
                                    [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                                    AR_TOF_Orbit_02(mzid)=STATS(1);                        
                                end                    
                            end
                            
                            if max(AR_TOF_Orbit_01)>=0.8 && max(AR_TOF_Orbit_02)>=0.8
                                ID3=ID3+1;
                                [Y01,I01]=max(AR_TOF_Orbit_01);
                                [Y02,I02]=max(AR_TOF_Orbit_02);
                                TOF_Orbit_goodalign0102_on_testing=[TOF_Orbit_goodalign0102_on_testing;i I01 I02];
                                Pepinformdetect(i).R_TOF_Orbit_Both_detect(k)=R_detect(k);
                                Pepinformdetect(i).C_TOF_Orbit_Both_detect(k)=C_detect(k);
                            else
                                Pepinformdetect(i).R_TOF_Orbit_Both_detect(k)=0;
                                Pepinformdetect(i).C_TOF_Orbit_Both_detect(k)=0;
                                if max(AR_TOF_Orbit_01)>=0.8
                                    ID1=ID1+1;
                                    [Y,I]=max(AR_TOF_Orbit_01);
                                    TOF_Orbit_goodalign01_on_testing=[TOF_Orbit_goodalign01_on_testing;i I];
                                    Pepinformdetect(i).R_TOF_Orbit_detect(k)=R_detect(k);
                                    Pepinformdetect(i).C_TOF_Orbit_detect(k)=C_detect(k);
                                else
                                    if max(AR_TOF_Orbit_02)>=0.8  
                                        ID2=ID2+1;
                                        [Y,I]=max(AR_TOF_Orbit_02);
                                        TOF_Orbit_goodalign02_on_testing=[TOF_Orbit_goodalign02_on_testing;i I];
                                        Pepinformdetect(i).R_TOF_Orbit_detect(k)=R_detect(k);
                                        Pepinformdetect(i).C_TOF_Orbit_detect(k)=C_detect(k);
                                    else
                                        Pepinformdetect(i).R_TOF_Orbit_detect(k)=0;
                                        Pepinformdetect(i).C_TOF_Orbit_detect(k)=0;
                                    end
                                end
                            end
                            clear TOF_XICs_matrix01 TOF_XICs_matrix02
                    end
                end
            end
            clear Orbit_XICs_matrix01 Orbit_XICs_matrix02
        else
            Pepinformdetect(i).ID_detect=0;
            Pepinformdetect(i).R_detect=0;
            Pepinformdetect(i).C_detect=0;
            Pepinformdetect(i).R_max_detect=0;
            Pepinformdetect(i).C_max_detect=0;            
        end
        
    else 
        Pepinformdetect(i).Orbit_align_JUD=0;
        AR_corr_Prob=betapdf(Testingdata{i}.AR,AR0102_PARMHAT1(1),AR0102_PARMHAT1(2));
        AR_non_corr_Prob=betapdf(Testingdata{i}.AR,AR0102_PARMHAT2(1),AR0102_PARMHAT2(2));
        AKL_corr_Prob=normpdf(Testingdata{i}.LOGAKL,AKL0102_MU_sig,AKL0102_SIGMA_sig);
        AKL_non_corr_Prob=normpdf(Testingdata{i}.LOGAKL,AKL0102_MU_notsig,AKL0102_SIGMA_notsig);
        AT_corr_Prob=normpdf(polyval(PP,Testingdata{i}.Normal_T01)-Testingdata{i}.Normal_T02,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
        AT_non_corr_Prob=normpdf(polyval(PP,Testingdata{i}.Normal_T01)-Testingdata{i}.Normal_T02,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
        Prob_AR_matrix=log(AR_corr_Prob./AR_non_corr_Prob);
        Prob_AKL_matrix=log(AKL_corr_Prob./AKL_non_corr_Prob);
        Prob_AT_matrix=log(AT_corr_Prob./AT_non_corr_Prob);
        Total_Prob_matrix=Prob_AR_matrix+Prob_AKL_matrix+Prob_AT_matrix;
        [R_detect,C_detect,V_detect]=find(Total_Prob_matrix>=0);
        [V_max_detect,I_max_detect]=max(Total_Prob_matrix);
        [V_max_detectv1,I_max_detectv1]=max(V_max_detect);    
        R_max_detect=I_max_detect(I_max_detectv1);
        C_max_detect=I_max_detectv1;
        if ~isempty(R_detect)
%             num=num+1;
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
end

length(unique(Orbit_goodalign_on_testing(:,1)))
length(unique(TOF_goodalign_corrBnon_on_testing))
length(unique(TOF_Orbit_goodalign0102_timescale_on_testing(:,1)))
length(unique(TOF_Orbit_goodalign01_on_testing(:,1)))
length(unique(TOF_Orbit_goodalign02_on_testing(:,1)))
length(unique(TOF_Orbit_goodalign0102_on_testing(:,1)))

save D:\Program\QTOF_replicate_identification\MATfile_v1\Testingdata Testingdata
save D:\Program\QTOF_replicate_identification\MATfile_v1\Pepinformdetect Pepinformdetect

%%%%%%%%%%%% if there are more than one peptide peak in the ms2 interval,
%%%%%%%%%%%% we need to choose one with highest R and T score
ID_TOF_Orbit_align0102intimescale=unique(TOF_Orbit_goodalign0102_timescale_on_testing(:,1));
num=0; ID_onepeak_inms2interval01=[];ID_onepeak_inms2interval02=[];
for i=1:length(Pepinformdetect)
    ID_F=find(ID_TOF_Orbit_align0102intimescale==i);
    if ~isempty(ID_F)
    
        KK=i;
        num=num+1;
        T_TOF_Orbit_Prob01=normpdf(polyval(PP_Normal_TOF2Orbit01,Testingdata{KK}.Normal_T_TOF01)-Testingdata{KK}.Normal_T_Orbit01,Mu_TOF2Orbit01,Sigma_TOF2Orbit01);
        R_TOF_Orbit_Prob01=gampdf(1-Testingdata{KK}.R_TOF01Orbit01,R_PHAT_TOF2Orbit01(1),R_PHAT_TOF2Orbit01(2));

        T_TOF_Orbit_Prob02=normpdf(polyval(PP_Normal_TOF2Orbit02,Testingdata{KK}.Normal_T_TOF02)-Testingdata{KK}.Normal_T_Orbit02,Mu_TOF2Orbit02,Sigma_TOF2Orbit02);
        R_TOF_Orbit_Prob02=gampdf(1-Testingdata{KK}.R_TOF02Orbit02,R_PHAT_TOF2Orbit02(1),R_PHAT_TOF2Orbit02(2));

        T_Detect_Prob01=T_TOF_Orbit_Prob01(Pepinformdetect(KK).R_TOF_Orbit_timescale_detect);
        R_Detect_Prob01=R_TOF_Orbit_Prob01(Pepinformdetect(KK).R_TOF_Orbit_timescale_detect);
        T_Detect_Prob02=T_TOF_Orbit_Prob02(Pepinformdetect(KK).C_TOF_Orbit_timescale_detect);
        R_Detect_Prob02=R_TOF_Orbit_Prob02(Pepinformdetect(KK).C_TOF_Orbit_timescale_detect);

        [Y,I]=max(log(T_Detect_Prob01)+log(R_Detect_Prob01)+log(T_Detect_Prob02)+log(R_Detect_Prob02));

        Pepinformdetect(KK).R_detect_basedon_timescale_RT=Pepinformdetect(KK).R_TOF_Orbit_timescale_detect(I);
        Pepinformdetect(KK).C_detect_basedon_timescale_RT=Pepinformdetect(KK).C_TOF_Orbit_timescale_detect(I);
        if length(T_Detect_Prob01)==1
            ID_onepeak_inms2interval01=[ID_onepeak_inms2interval01; i];            
        end
        if length(T_Detect_Prob02)==1
            ID_onepeak_inms2interval02=[ID_onepeak_inms2interval02; i];            
        end
    else 
        KK=i;
        Pepinformdetect(KK).R_detect_basedon_timescale_RT=0;
        Pepinformdetect(KK).C_detect_basedon_timescale_RT=0;
  
        
    end
end
ID_onepeak_inms2interval0102=intersect(ID_onepeak_inms2interval01,ID_onepeak_inms2interval02);
%%%%%%%%%%%%

%%%%%%%%%%%% generate Align Result matrix
load XUEPO_Orbi_5569_f4_info_newtxt
for i=2:size(matrix,1)
    ID_infinal=[str2num(matrix{i,9}), str2num(matrix{i,17}), str2num(matrix{i,25})];
    ID_data_index=find(ID_infinal~=0);
    switch ID_data_index(1)
        case 1
            Protein_matrix_final{i-1}= protein01_final{ID_infinal(ID_data_index(1))};           
        case 2
            Protein_matrix_final{i-1}= protein02_final{ID_infinal(ID_data_index(1))}; 
        case 3
            Protein_matrix_final{i-1}= protein03_final{ID_infinal(ID_data_index(1))}; 
    end
end
save D:\Program\QTOF_replicate_identification\MATfile_v1\Protein_matrix_final Protein_matrix_final
matrix_new=matrix(Com_ID_0102+1,:);
Protein_matrix_final_new=Protein_matrix_final(Com_ID_0102);
Testing_align_id=unique(TOF_Orbit_goodalign0102_timescale_on_testing(:,1));

matrix_TOF{1,1}='pepseq';
matrix_TOF{1,2}='cs01';matrix_TOF{1,3}='mass01';matrix_TOF{1,4}='timepoint01';matrix_TOF{1,5}='remarks01';
matrix_TOF{1,6}='INTstart01';matrix_TOF{1,7}='INTend01';matrix_TOF{1,8}='Prob01';matrix_TOF{1,9}='INDindata01';
matrix_TOF{1,10}='cs02';matrix_TOF{1,11}='mass02';matrix_TOF{1,12}='timepoint02';matrix_TOF{1,13}='remarks02';
matrix_TOF{1,14}='INTstart02';matrix_TOF{1,15}='INTend02';matrix_TOF{1,16}='Prob02';matrix_TOF{1,17}='INDindata02';
for i=1:length(Training_id)
    
    iso_matrix_TOF(i,:)=iso_com(Training_id(i),:);
    Protein_matrix_TOF{i}=Protein_matrix_final_new{Training_id(i)};
    
    matrix_TOF{i+1,1}=matrix_new{Training_id(i),1};
    matrix_TOF{i+1,2}=matrix_new{Training_id(i),2};
    matrix_TOF{i+1,3}=matrix_new{Training_id(i),3};
    matrix_TOF{i+1,4}=0;
    matrix_TOF{i+1,5}=1;
    Training_Corr_peak_id01=Training_Row_Column_Num(i,1);
    matrix_TOF{i+1,6}=IntervalList_com01{Training_id(i)}.intervallist_after7_combine(Training_Corr_peak_id01,1);
    matrix_TOF{i+1,7}=IntervalList_com01{Training_id(i)}.intervallist_after7_combine(Training_Corr_peak_id01,2);
    matrix_TOF{i+1,8}=0;
    matrix_TOF{i+1,9}=Training_id(i);
    matrix_TOF{i+1,10}=matrix_new{Training_id(i),10};
    matrix_TOF{i+1,11}=matrix_new{Training_id(i),11};
    matrix_TOF{i+1,12}=0;
    matrix_TOF{i+1,13}=1;
    Training_Corr_peak_id02=Training_Row_Column_Num(i,2);
    matrix_TOF{i+1,14}=IntervalList_com02{Training_id(i)}.intervallist_after7_combine(Training_Corr_peak_id02,1);
    matrix_TOF{i+1,15}=IntervalList_com02{Training_id(i)}.intervallist_after7_combine(Training_Corr_peak_id02,2);
    matrix_TOF{i+1,16}=0;
    matrix_TOF{i+1,17}=Training_id(i);
end
for i=1:length(Testing_align_id) 

    mm=Testing_align_id(i);
    
    iso_matrix_TOF(i+length(Training_id),:)=iso_com(Testing_id(mm),:);
    Protein_matrix_TOF{i+length(Training_id)}=Protein_matrix_final_new{Testing_id(mm)};
    
    matrix_TOF{i+1+length(Training_id),1}=matrix_new{Testing_id(mm),1};
    matrix_TOF{i+1+length(Training_id),2}=matrix_new{Testing_id(mm),2};
    matrix_TOF{i+1+length(Training_id),3}=matrix_new{Testing_id(mm),3};
    matrix_TOF{i+1+length(Training_id),4}=0;
    matrix_TOF{i+1+length(Training_id),5}=1;
    Testing_Corr_peak_id01=Pepinformdetect(mm).R_detect_basedon_timescale_RT;
    matrix_TOF{i+1+length(Training_id),6}=IntervalList_com01{Testing_id(mm)}.intervallist_after7_combine(Testing_Corr_peak_id01,1);
    matrix_TOF{i+1+length(Training_id),7}=IntervalList_com01{Testing_id(mm)}.intervallist_after7_combine(Testing_Corr_peak_id01,2);
    matrix_TOF{i+1+length(Training_id),8}=0;
    matrix_TOF{i+1+length(Training_id),9}=Testing_id(mm);
    matrix_TOF{i+1+length(Training_id),10}=matrix_new{Testing_id(mm),10};
    matrix_TOF{i+1+length(Training_id),11}=matrix_new{Testing_id(mm),11};
    matrix_TOF{i+1+length(Training_id),12}=0;
    matrix_TOF{i+1+length(Training_id),13}=1;
    Testing_Corr_peak_id02=Pepinformdetect(mm).C_detect_basedon_timescale_RT;
    matrix_TOF{i+1+length(Training_id),14}=IntervalList_com02{Testing_id(mm)}.intervallist_after7_combine(Testing_Corr_peak_id02,1);
    matrix_TOF{i+1+length(Training_id),15}=IntervalList_com02{Testing_id(mm)}.intervallist_after7_combine(Testing_Corr_peak_id02,2);
    matrix_TOF{i+1+length(Training_id),16}=0;
    matrix_TOF{i+1+length(Training_id),17}=Testing_id(mm);
end

size(iso_matrix_TOF)
size(Protein_matrix_TOF)
save D:\Program\QTOF_replicate_identification\MATfile_v1\TOF5574_f4_0102_information matrix_TOF iso_matrix_TOF Protein_matrix_TOF
save D:\Program\QTOF_replicate_identification\MATfile_v1\TOF5574_f4_0102_XICs Orbit_XICs_com01 TOF_XICs_com01 Orbit_XICs_com02 TOF_XICs_com02
%%%%%%%%%%%%
%%%%%%%%%% Check mass diff

Training_Matrix_TOF=matrix_TOF(2:length(Training_id)+1,:);
Testing_Matrix_TOF=matrix_TOF(length(Training_id)+1:end,:);
for i=1:size(Training_Matrix_TOF,1)
    PEP_spec=Training_Matrix_TOF{i,1};
    cs01=str2num(Training_Matrix_TOF{i,2});
    mass01=str2num(Training_Matrix_TOF{i,3});
    mz01=(mass01+cs01*1.0073)/cs01;
    Interval_start01=Training_Matrix_TOF{i,6};
    Interval_end01=Training_Matrix_TOF{i,7};
    PPM=20;
    MZ_interval=[mz01-mz01*PPM/10^6,mz01+mz01*PPM/10^6];
    [Int_Matrix_mz_time01,mz_unique01,Pep_mz]=Generate_Pep_Mz_Time_matrix(TOF_peakl01,MZ_interval,Interval_start01,Interval_end01);
    X_index=mz_unique01*ones(1,size(Int_Matrix_mz_time01,2));
    Y_index=ones(size(Int_Matrix_mz_time01,1),1)*[Interval_start01:Interval_end01];
    Z_index=Int_Matrix_mz_time01;
    figure
    stem3(X_index,Y_index,Z_index)
    grid on
end
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% check the detected intervals
%%%%%%%%%% plot Orbit (01 & 02) and TOF (01 & 02) XICs

TOF_Orbit_on_testing=unique(TOF_Orbit_goodalign0102_timescale_on_testing(:,1));
TOF_Orbit_notalign_on_testing=setdiff(Orbit_goodalign_on_testing(:,1),TOF_Orbit_on_testing);
% TOF_Orbit_on_testing=unique(TOF_Orbit_goodalign01_on_testing(:,1));
% TOF_Orbit_notalign_on_testing=setdiff(unique(TOF_Orbit_goodalign0102_timescale_on_testing(:,1)),TOF_Orbit_on_testing);
% TOF_Orbit_on_testing=unique(TOF_Orbit_goodalign02_on_testing(:,1));
% TOF_Orbit_notalign_on_testing=setdiff(unique(TOF_Orbit_goodalign0102_timescale_on_testing(:,1)),TOF_Orbit_on_testing);
% TOF_Orbit_on_testing=unique(TOF_Orbit_goodalign0102_on_testing(:,1));
% TOF_Orbit_notalign_on_testing=setdiff(unique(TOF_Orbit_goodalign0102_timescale_on_testing(:,1)),TOF_Orbit_on_testing);

for kkk=601:605%length(Pepinformdetect)
%     i=TOF_Orbit_notalign_on_testing(kkk);
    i=TOF_Orbit_on_testing(kkk);
    ID_detect=Pepinformdetect(i).ID_detect;
%     R_detect=Pepinformdetect(i).R_detect;
%     C_detect=Pepinformdetect(i).C_detect;
%     R_detect=Pepinformdetect(i).R_max_detect;
%     C_detect=Pepinformdetect(i).C_max_detect;
%     R_detect=Pepinformdetect(i).R_TOF_Orbit_timescale_detect;
%     C_detect=Pepinformdetect(i).C_TOF_Orbit_timescale_detect;
%     R_detect=Pepinformdetect(i).R_TOF_Orbit_detect;
%     C_detect=Pepinformdetect(i).C_TOF_Orbit_detect;
%     R_detect=Pepinformdetect(i).R_TOF_Orbit_Both_detect;
%     C_detect=Pepinformdetect(i).C_TOF_Orbit_Both_detect;

    R_detect_basedon_timescale_RT=Pepinformdetect(i).R_detect_basedon_timescale_RT;
    C_detect_basedon_timescale_RT=Pepinformdetect(i).C_detect_basedon_timescale_RT;
    
    R_detectv1=R_detect_basedon_timescale_RT(R_detect_basedon_timescale_RT+C_detect_basedon_timescale_RT~=0);
    C_detectv1=C_detect_basedon_timescale_RT(R_detect_basedon_timescale_RT+C_detect_basedon_timescale_RT~=0);
    
    R_detect=Pepinformdetect(i).R_detect;
    C_detect=Pepinformdetect(i).C_detect;
    
    R_intimescale_detect=Pepinformdetect(i).R_TOF_Orbit_timescale_detect;
    C_intimescale_detect=Pepinformdetect(i).C_TOF_Orbit_timescale_detect;
    
    Interval_index_M01=TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_intimescale_detect,:))/max(TOF_retentiont01l1);
    Interval_index_M02=TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_intimescale_detect,:))/max(TOF_retentiont02l1);
    
    

    
    %%%%%
    Orbit_Time01=Testing_Orbit_AMTms2time_com(ID_detect,1);    
    MS2_scan_number01=find(retentiont01l1>=Orbit_Time01);
    Orbit_detect_scan01=Testing_Orbit_AMTms2time_com(ID_detect,2:3);
    Orbit_Time02=Testing_Orbit_AMTms2time_com(ID_detect,4);    
    MS2_scan_number02=find(retentiont02l1>=Orbit_Time02);
    Orbit_detect_scan02=Testing_Orbit_AMTms2time_com(ID_detect,5:6); 
    %%%%%
    TOF_instrument_h_int=10^7;    TOF_de_range=10^5;    TOF_Times_noise_std=3;
    Orbit_instrument_h_int=10^7;    Orbit_de_range=10^4;    Orbit_Times_noise_std=6;
    for j=1:6
        TOF_XICs01_matrix(:,j)=TOF_XICs_testing01{j}(:,ID_detect);
        TOF_Threshold01(j)=getThreshold_range(TOF_XICs01_matrix(:,j),TOF_instrument_h_int,TOF_de_range,TOF_Times_noise_std);
        TOF_XICs02_matrix(:,j)=TOF_XICs_testing02{j}(:,ID_detect);
        TOF_Threshold02(j)=getThreshold_range(TOF_XICs02_matrix(:,j),TOF_instrument_h_int,TOF_de_range,TOF_Times_noise_std);
        Orbit_XICs01_matrix(:,j)=Orbit_XICs_testing01{j}(:,ID_detect);
        Orbit_Threshold01(j)=getThreshold_range(Orbit_XICs01_matrix(:,j),Orbit_instrument_h_int,Orbit_de_range,Orbit_Times_noise_std);
        Orbit_XICs02_matrix(:,j)=Orbit_XICs_testing02{j}(:,ID_detect);
        Orbit_Threshold02(j)=getThreshold_range(Orbit_XICs02_matrix(:,j),Orbit_instrument_h_int,Orbit_de_range,Orbit_Times_noise_std);
       
    end
 
    colorarray=['r', 'k', 'g', 'b', 'm', 'y'];

    TOF_time_max01=max(TOF_retentiont01l1);
    Orbit_time_max01=max(retentiont01l1);
    
    figure
    subplot(2,2,1)
    for j=1:6
        plot(TOF_retentiont01l1/TOF_time_max01,TOF_XICs01_matrix(:,j),colorarray(j))
        hold on
        plot(TOF_retentiont01l1/TOF_time_max01,TOF_Threshold01(j)*ones(1,length(TOF_retentiont01l1)),colorarray(j))
    end
    text(0,max(max(TOF_XICs01_matrix)),num2str(Interval_index_M01))
    height=1000;
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detectv1,1))/TOF_time_max01,2*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detectv1,1)),1),'r*')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detectv1,2))/TOF_time_max01,2*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detectv1,2)),1),'k*')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_intimescale_detect,1))/TOF_time_max01,1.75*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_intimescale_detect,1)),1),'rv')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_intimescale_detect,2))/TOF_time_max01,1.75*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_intimescale_detect,2)),1),'kv')        
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detect,1))/TOF_time_max01,1.5*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detect,1)),1),'rh')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detect,2))/TOF_time_max01,1.5*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detect,2)),1),'kh')        
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(:,1))/TOF_time_max01,1.25*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(:,1)),1),'ro')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine(:,2))/TOF_time_max01,1.25*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine(:,2)),1),'ko')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine_beforecheckL(:,1))/TOF_time_max01,1*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine_beforecheckL(:,1)),1),'rd')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7_combine_beforecheckL(:,2))/TOF_time_max01,1*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7_combine_beforecheckL(:,2)),1),'kd')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7(:,1))/TOF_time_max01,0.75*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7(:,1)),1),'rx')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after7(:,2))/TOF_time_max01,0.75*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after7(:,2)),1),'kx')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after6(:,1))/TOF_time_max01,0.50*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after6(:,1)),1),'rs')
    stem(TOF_retentiont01l1(IntervalList_testing01{ID_detect}.intervallist_after6(:,2))/TOF_time_max01,0.50*height*ones(length(IntervalList_testing01{ID_detect}.intervallist_after6(:,2)),1),'ks')
    grid on
    
    subplot(2,2,2)
    for j=1:6
        plot(retentiont01l1/Orbit_time_max01,Orbit_XICs01_matrix(:,j),colorarray(j))
        hold on
        plot(retentiont01l1/Orbit_time_max01,Orbit_Threshold01(j)*ones(1,length(retentiont01l1)),colorarray(j))
    end
    text(0,max(max(Orbit_XICs01_matrix)),num2str(retentiont01l1(Orbit_detect_scan01)/Orbit_time_max01))
    stem(retentiont01l1(MS2_scan_number01(1))/Orbit_time_max01,20000000*ones(length(MS2_scan_number01(1)),1),'g*')
    stem(retentiont01l1(Orbit_detect_scan01(1))/Orbit_time_max01,10000000*ones(length(Orbit_detect_scan01(1)),1),'ro')
    stem(retentiont01l1(Orbit_detect_scan01(2))/Orbit_time_max01,10000000*ones(length(Orbit_detect_scan01(2)),1),'ko')
    grid on
    
    TOF_time_max02=max(TOF_retentiont02l1);
    Orbit_time_max02=max(retentiont02l1);
    
    subplot(2,2,3)
    for j=1:6
        plot(TOF_retentiont02l1/TOF_time_max02,TOF_XICs02_matrix(:,j),colorarray(j))
        hold on
        plot(TOF_retentiont02l1/TOF_time_max02,TOF_Threshold02(j)*ones(1,length(TOF_retentiont02l1)),colorarray(j))
    end
    text(0,max(max(TOF_XICs02_matrix)),num2str(Interval_index_M02))
    height=1000;
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detectv1,1))/TOF_time_max02,2*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detectv1,1)),1),'r*')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detectv1,2))/TOF_time_max02,2*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detectv1,2)),1),'k*')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_intimescale_detect,1))/TOF_time_max01,1.75*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_intimescale_detect,1)),1),'rv')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_intimescale_detect,2))/TOF_time_max01,1.75*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_intimescale_detect,2)),1),'kv')        
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detect,1))/TOF_time_max01,1.5*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detect,1)),1),'rh')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detect,2))/TOF_time_max01,1.5*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detect,2)),1),'kh')        
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(:,1))/TOF_time_max02,1.25*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(:,1)),1),'ro')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine(:,2))/TOF_time_max02,1.25*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine(:,2)),1),'ko')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine_beforecheckL(:,1))/TOF_time_max02,1*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine_beforecheckL(:,1)),1),'rd')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7_combine_beforecheckL(:,2))/TOF_time_max02,1*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7_combine_beforecheckL(:,2)),1),'kd')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7(:,1))/TOF_time_max02,0.75*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7(:,1)),1),'rx')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after7(:,2))/TOF_time_max02,0.75*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after7(:,2)),1),'kx')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after6(:,1))/TOF_time_max02,0.50*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after6(:,1)),1),'rs')
    stem(TOF_retentiont02l1(IntervalList_testing02{ID_detect}.intervallist_after6(:,2))/TOF_time_max02,0.50*height*ones(length(IntervalList_testing02{ID_detect}.intervallist_after6(:,2)),1),'ks')
    grid on
    
    subplot(2,2,4)
    for j=1:6
        plot(retentiont02l1/Orbit_time_max02,Orbit_XICs02_matrix(:,j),colorarray(j))
        hold on
        plot(retentiont02l1/Orbit_time_max02,Orbit_Threshold02(j)*ones(1,length(retentiont02l1)),colorarray(j))
    end
    text(0,max(max(Orbit_XICs02_matrix)),num2str(retentiont02l1(Orbit_detect_scan02)/Orbit_time_max02))
    stem(retentiont02l1(MS2_scan_number02(1))/Orbit_time_max02,20000000*ones(length(MS2_scan_number02(1)),1),'g*')
    stem(retentiont02l1(Orbit_detect_scan02(1))/Orbit_time_max02,10000000*ones(length(Orbit_detect_scan02(1)),1),'ro')
    stem(retentiont02l1(Orbit_detect_scan02(2))/Orbit_time_max02,10000000*ones(length(Orbit_detect_scan02(2)),1),'ko')
    grid on
    

    %%%%%%%%%%%%%%%%%%%% TOF spec pep information
    index=Testing_id(ID_detect);
    
    TOF_spec_mass=TOF_mass_com(index);
    TOF_spec_cs=TOF_charge_com(index);
    TOF_spec_time=TOF_ms2information_com(index);
    mzvalue=(TOF_spec_mass+TOF_spec_cs*1.0073)/TOF_spec_cs;
    mzlist=[mzvalue mzvalue+1.00335/csvalue mzvalue+2.00547/csvalue mzvalue+3.00882/csvalue mzvalue+4.00849/csvalue mzvalue+5.01184/csvalue];
    %%%%%%%%%%%%%%%%%%%%
    
    Interval_matrix01=IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detect,:);
    MZ_mono_interval(1,:)=mzlist-mzlist*80/10^6;
    MZ_mono_interval(2,:)=mzlist+mzlist*80/10^6;
    for j=1:size(Interval_matrix01,1)
        [Int_Matrix_mz_time_mono,mz_unique_mono,Pep_mz_mono]=Generate_Pep_Mz_Time_matrix(TOF_peakl01,MZ_mono_interval(:,1),Interval_matrix01(j,1),Interval_matrix01(j,2));
        [Int_Matrix_mz_time_iso1st,mz_unique_iso1st,Pep_mz_iso1st]=Generate_Pep_Mz_Time_matrix(TOF_peakl01,MZ_mono_interval(:,2),Interval_matrix01(j,1),Interval_matrix01(j,2));
        
        xxx=mz_unique_mono*ones(1,size(Int_Matrix_mz_time_mono,2));
        yyy=ones(size(Int_Matrix_mz_time_mono,1),1)*[Interval_matrix01(j,1):Interval_matrix01(j,2)];
        xxx1=mzlist(1)*ones(1,size(Int_Matrix_mz_time_mono,2));
        yyy1=[Interval_matrix01(j,1):Interval_matrix01(j,2)];
        zzz1=zeros(1,size(Int_Matrix_mz_time_mono,2));
        figure
%         surf(xxx,yyy,Int_Matrix_mz_time_mono)
%         hold on
    
        stem3(xxx,TOF_retentiont01l1(yyy)/TOF_time_max01,Int_Matrix_mz_time_mono,'k.')
        hold on
        plot3(xxx1,TOF_retentiont01l1(yyy1)/TOF_time_max01,zzz1,'r');
        grid on
        X_text=double(round(min(min(xxx))*100)/100);
        Y_text=double(round(min(min(TOF_retentiont01l1(yyy)/TOF_time_max01))*100)/100);
        text(X_text,Y_text,0,num2str(kkk));
        
%         xxx=mz_unique_mono*ones(1,size(Int_Matrix_mz_time_iso1st,2));
%         yyy=ones(size(Int_Matrix_mz_time_iso1st,1),1)*[Interval_matrix01(j,1):Interval_matrix01(j,2)];
%         figure
%         stem3(xxx,yyy,Int_Matrix_mz_time_iso1st)
        
      
    end
    Mass_standard=TOF_spec_mass;  


    
    
    
end
%%%%%%%%%%



% %%%%%%%%%% check spec peptide
% i=101;
% %%%%%%%%%% check spec pep TOF detected intervals 
% ID_detect=Pepinformdetect(i).ID_detect;
% R_detect=Pepinformdetect(i).R_detect;
% C_detect=Pepinformdetect(i).C_detect;
% %%%%%%%%%% check the interval list
% % IntervalList_testing01{ID_detect}.intervallist_after7
% % IntervalList_testing02{ID_detect}.intervallist_after7
% % IntervalList_testing01{ID_detect}.intervallist_after7_combine_beforecheckL
% % IntervalList_testing02{ID_detect}.intervallist_after7_combine_beforecheckL
% % IntervalList_testing01{ID_detect}.intervallist_after7_combine
% % IntervalList_testing02{ID_detect}.intervallist_after7_combine
% % IntervalList_testing01{ID_detect}.intervallist_after7_combine(R_detect,:)
% % IntervalList_testing02{ID_detect}.intervallist_after7_combine(C_detect,:)
% %%%%%%%%%%
% %%%%%%%%%% check XICs
% Spec_pep=Pep_com(Testing_id(ID_detect));
% 
% [peptidenew, massdiffList, isHeavy]=modprocess({Spec_pep{1}});
% peptidenew=peptidenew{1};
% [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(Spec_pep{1}));
% ID_cs=find(Orbit_AMTcs_com(Testing_id(ID_detect),:)>0);
% csvalue=Orbit_AMTcs_com(Testing_id(ID_detect),ID_cs(1));
% mzvalue=(weight+csvalue*1.0073)/csvalue;
% Spec_mzList=[mzvalue mzvalue+1.00335/csvalue mzvalue+2.00547/csvalue mzvalue+3.00882/csvalue mzvalue+4.00849/csvalue mzvalue+5.01184/csvalue];
% 
% load ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_01_f4.mzXML.centroid1.peak.mat
% peakl01_cent1=peakl;
% clear peakl
% 
% for tolerance=20%:10:100;
%     for j=1:6
% %         TOF_spec_XICs_total(:,j)=getXICs(peakl01_cent1,Spec_mzList(j),tolerance);
%         TOF_spec_XICs_centroid0_total(:,j)=getXICs(peakl01,Spec_mzList(j),tolerance);
%     end
%     for j=1:6
% %         TOF_spec_XICs(:,j)=TOF_spec_XICs_total(1:length(TOF_retentiont01l1),j);
%         TOF_spec_XICs_centroid0(:,j)=TOF_spec_XICs_centroid0_total(1:length(TOF_retentiont01l1),j);
%         instrument_h_int=10^7;
%         de_range=10^6;
%         Times_noise_std=3;
% %         Interval_int_threshold(j)=TOF_getThreshold_range(TOF_spec_XICs(:,j),instrument_h_int, de_range,Times_noise_std);
%         Interval_int_threshold_centroid0(j)=TOF_getThreshold_range(TOF_spec_XICs_centroid0(:,j),instrument_h_int, de_range,Times_noise_std);
%     end
%     
%     %%%%%%%%%%%%% interval detection on one peptide
%     iso_spec_pep=iso_testing(ID_detect,:);
%     if Testing_Orbit_AMTms2time_com(ID_detect,1)~=0
%         Orbit_MS2_time=Testing_Orbit_AMTms2time_com(ID_detect,1);
%     else
%         Orbit_MS2_time=sum(retentiont01l1(Testing_Orbit_AMTms2time_com(ID_detect,2:3)))/2;
%     end
%     [IntervalList_cent0,Jud_detect_good_cent0]=TOF_Verification_intervaldetection_Spec_Pep(Spec_pep,iso_spec_pep,Orbit_MS2_time,TOF_spec_XICs_centroid0,TOF_retentiont01l1);
% %     [IntervalList_cent1,Jud_detect_good_cent1]=TOF_Verification_intervaldetection_Spec_Pep(Spec_pep,iso_spec_pep,Orbit_MS2_time,TOF_spec_XICs,TOF_retentiont01l1);
%     
%     %%%%%%%%%%%%%
%     
%     colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
% %     figure
% %     for j=1:6
% %         plot(TOF_retentiont01l1,TOF_spec_XICs(:,j),colorarray(j))
% %         hold on
% %         plot(TOF_retentiont01l1,Interval_int_threshold(j)*ones(1,length(TOF_retentiont01l1)),colorarray(j))
% %     end
% %     height=1000;
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after7_combine(:,1)),1.25*height*ones(length(IntervalList_cent1.intervallist_after7_combine(:,1)),1),'ro')
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after7_combine(:,2)),1.25*height*ones(length(IntervalList_cent1.intervallist_after7_combine(:,2)),1),'ko')
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after7_combine_beforecheckL(:,1)),1*height*ones(length(IntervalList_cent1.intervallist_after7_combine_beforecheckL(:,1)),1),'rd')
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after7_combine_beforecheckL(:,2)),1*height*ones(length(IntervalList_cent1.intervallist_after7_combine_beforecheckL(:,2)),1),'kd')
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after7(:,1)),0.75*height*ones(length(IntervalList_cent1.intervallist_after7(:,1)),1),'rx')
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after7(:,2)),0.75*height*ones(length(IntervalList_cent1.intervallist_after7(:,2)),1),'kx')
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after6(:,1)),0.50*height*ones(length(IntervalList_cent1.intervallist_after6(:,1)),1),'rs')
% %     stem(TOF_retentiont01l1(IntervalList_cent1.intervallist_after6(:,2)),0.50*height*ones(length(IntervalList_cent1.intervallist_after6(:,2)),1),'ks')
% %     grid on;
% 
%     figure
%     for j=1:6
%         plot(TOF_retentiont01l1,TOF_spec_XICs_centroid0(:,j),colorarray(j))
%         hold on
%         plot(TOF_retentiont01l1,Interval_int_threshold_centroid0(j)*ones(1,length(TOF_retentiont01l1)),colorarray(j))
%     end
%     height=1000;
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after7_combine(:,1)),1.25*height*ones(length(IntervalList_cent0.intervallist_after7_combine(:,1)),1),'ro')
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after7_combine(:,2)),1.25*height*ones(length(IntervalList_cent0.intervallist_after7_combine(:,2)),1),'ko')
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after7_combine_beforecheckL(:,1)),1*height*ones(length(IntervalList_cent0.intervallist_after7_combine_beforecheckL(:,1)),1),'rd')
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after7_combine_beforecheckL(:,2)),1*height*ones(length(IntervalList_cent0.intervallist_after7_combine_beforecheckL(:,2)),1),'kd')
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after7(:,1)),0.75*height*ones(length(IntervalList_cent0.intervallist_after7(:,1)),1),'rx')
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after7(:,2)),0.75*height*ones(length(IntervalList_cent0.intervallist_after7(:,2)),1),'kx')
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after6(:,1)),0.50*height*ones(length(IntervalList_cent0.intervallist_after6(:,1)),1),'rs')
%     stem(TOF_retentiont01l1(IntervalList_cent0.intervallist_after6(:,2)),0.50*height*ones(length(IntervalList_cent0.intervallist_after6(:,2)),1),'ks')
%     grid on
% 
% end
% %%%%%%%%%%



%%%%%%%%%%%%%%%



