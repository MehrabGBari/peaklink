clc
clear all

%%%%%%%%%%%%%%%%% load Orbitrap ms2 information
%%%%%%%% XICs01 02 03 totalmzList finalchargeList matrix and so on
ii=5569;
jj=4;
eval(['load ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_01_f', num2str(jj), '.mzXML.centroid0.peak.mat'])
Orbi_data01l1=peakl;
Orbi_retentiont01_l1l2=retentiont;
Orbi_id_level1_01=find(Orbi_retentiont01_l1l2(:,2)==1);
Orbi_retentiont01l1=Orbi_retentiont01_l1l2(Orbi_id_level1_01,1);
clear peakl retentiont

eval(['load ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_02_f', num2str(jj), '.mzXML.centroid0.peak.mat'])
Orbi_data02l1=peakl;
Orbi_retentiont02_l1l2=retentiont;
Orbi_id_level1_02=find(Orbi_retentiont02_l1l2(:,2)==1);
Orbi_retentiont02l1=Orbi_retentiont02_l1l2(Orbi_id_level1_02,1);
clear peakl retentiont

eval(['load ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_03_f', num2str(jj), '.mzXML.centroid0.peak.mat'])
Orbi_data03l1=peakl;
Orbi_retentiont03_l1l2=retentiont;
Orbi_id_level1_03=find(Orbi_retentiont03_l1l2(:,2)==1);
Orbi_retentiont03l1=Orbi_retentiont03_l1l2(Orbi_id_level1_03,1);
clear peakl retentiont

% file_xls1=['ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_01_f',num2str(jj),'.tandem.pep.xls'];
% file_xls2=['ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_02_f',num2str(jj),'.tandem.pep.xls'];
% file_xls3=['ZHA_27_',num2str(ii),'_21APR11_CELL_VEL_HUM_TT_JL_2D_03_f',num2str(jj),'.tandem.pep.xls'];
% 
% Orbi_data01_information=Saveinformation_1_6data(file_xls1);%% read in xls for 01data
% Orbi_data02_information=Saveinformation_1_6data(file_xls2);%% read in xls for 02data
% Orbi_data03_information=Saveinformation_1_6data(file_xls3);%% read in xls for 03data
% 
% save Orbi_data01_information Orbi_data01_information
% save Orbi_data02_information Orbi_data02_information
% save Orbi_data03_information Orbi_data03_information

load Orbi_data01_information
load Orbi_data02_information
load Orbi_data03_information

Orbi_data01_information.remarks=1*ones(length(Orbi_data01_information.peptide),1);
Orbi_data02_information.remarks=2*ones(length(Orbi_data02_information.peptide),1);
Orbi_data03_information.remarks=3*ones(length(Orbi_data03_information.peptide),1);

[Orbi_data01v1,I01,P01_once]=filteroverlapp(Orbi_data01_information);
[Orbi_data02v1,I02,P02_once]=filteroverlapp(Orbi_data02_information);
[Orbi_data03v1,I03,P03_once]=filteroverlapp(Orbi_data03_information);

ID_XIC01=I01(P01_once);
ID_XIC02=I02(P02_once);
ID_XIC03=I03(P03_once);

P_threshold=0.9; %%% set up the threshold of Prophet Score
Id01_95=find(Orbi_data01v1.probability>=P_threshold);
Id02_95=find(Orbi_data02v1.probability>=P_threshold);
Id03_95=find(Orbi_data03v1.probability>=P_threshold);

%%%%%%%%% get basic information of data01 data02 data03
Orbi_pep01_original=Orbi_data01v1.peptide;
Orbi_chargestate01_original=Orbi_data01v1.assumed_charge;
Orbi_ms2information01_original=Orbi_data01v1.retention_time_sec;
Orbi_mass01_original=Orbi_data01v1.calc_neutral_pep_mass;
Orbi_remarks01_original=Orbi_data01v1.remarks;
Orbi_Probability01_original=Orbi_data01v1.probability;
Orbi_mzratio01_original=Orbi_data01v1.MZratio;
Orbi_protein01_original=Orbi_data01v1.protein;

Orbi_pep02_original=Orbi_data02v1.peptide;
Orbi_chargestate02_original=Orbi_data02v1.assumed_charge;
Orbi_ms2information02_original=Orbi_data02v1.retention_time_sec;
Orbi_mass02_original=Orbi_data02v1.calc_neutral_pep_mass;
Orbi_remarks02_original=Orbi_data02v1.remarks;
Orbi_Probability02_original=Orbi_data02v1.probability;
Orbi_mzratio02_original=Orbi_data02v1.MZratio;
Orbi_protein02_original=Orbi_data02v1.protein;

Orbi_pep03_original=Orbi_data03v1.peptide;
Orbi_chargestate03_original=Orbi_data03v1.assumed_charge;
Orbi_ms2information03_original=Orbi_data03v1.retention_time_sec;
Orbi_mass03_original=Orbi_data03v1.calc_neutral_pep_mass;
Orbi_remarks03_original=Orbi_data03v1.remarks;
Orbi_Probability03_original=Orbi_data03v1.probability;
Orbi_mzratio03_original=Orbi_data03v1.MZratio;
Orbi_protein03_original=Orbi_data03v1.protein;
%%%%%%%%%%

%%%%%%%%%% get pep information that has Prophet Score bigger than 0.95
Orbi_pep01=Orbi_pep01_original(Id01_95);
Orbi_chargestate01=Orbi_chargestate01_original(Id01_95);
Orbi_ms2information01=Orbi_ms2information01_original(Id01_95);
Orbi_mass01=Orbi_mass01_original(Id01_95);
Orbi_remarks01=Orbi_remarks01_original(Id01_95);
Orbi_Probability01=Orbi_Probability01_original(Id01_95);
Orbi_protein01=Orbi_protein01_original(Id01_95);

Orbi_pep02=Orbi_pep02_original(Id02_95);
Orbi_chargestate02=Orbi_chargestate02_original(Id02_95);
Orbi_ms2information02=Orbi_ms2information02_original(Id02_95);
Orbi_mass02=Orbi_mass02_original(Id02_95);
Orbi_remarks02=Orbi_remarks02_original(Id02_95);
Orbi_Probability02=Orbi_Probability02_original(Id02_95);
Orbi_protein02=Orbi_protein02_original(Id02_95);

Orbi_pep03=Orbi_pep03_original(Id03_95);
Orbi_chargestate03=Orbi_chargestate03_original(Id03_95);
Orbi_ms2information03=Orbi_ms2information03_original(Id03_95);
Orbi_mass03=Orbi_mass03_original(Id03_95);
Orbi_remarks03=Orbi_remarks03_original(Id03_95);
Orbi_Probability03=Orbi_Probability03_original(Id03_95);
Orbi_protein03=Orbi_protein03_original(Id03_95);

%%%%%%%%%%%%%%%%%%%%%% only calculate the Orbitrap 01 02 overlapped
%%%%%%%%%%%%%%%%%%%%%% peptides
[Pepcommon0102, IA01, IB02]=intersect(Orbi_pep01,Orbi_pep02);
ID_same_cs=find(Orbi_chargestate01(IA01)-Orbi_chargestate02(IB02)==0);

Pepcommon0102v1=Pepcommon0102(ID_same_cs);

Orbi_pep01v1=Orbi_pep01(IA01(ID_same_cs));
Orbi_chargestate01v1=Orbi_chargestate01(IA01(ID_same_cs));
Orbi_ms2information01v1=Orbi_ms2information01(IA01(ID_same_cs));
Orbi_mass01v1=Orbi_mass01(IA01(ID_same_cs));
Orbi_remarks01v1=Orbi_remarks01(IA01(ID_same_cs));
Orbi_Probability01v1=Orbi_Probability01(IA01(ID_same_cs));
Orbi_protein01v1=Orbi_protein01(IA01(ID_same_cs));

Orbi_pep02v1=Orbi_pep02(IB02(ID_same_cs));
Orbi_chargestate02v1=Orbi_chargestate02(IB02(ID_same_cs));
Orbi_ms2information02v1=Orbi_ms2information02(IB02(ID_same_cs));
Orbi_mass02v1=Orbi_mass02(IB02(ID_same_cs));
Orbi_remarks02v1=Orbi_remarks02(IB02(ID_same_cs));
Orbi_Probability02v1=Orbi_Probability02(IB02(ID_same_cs));
Orbi_protein02v1=Orbi_protein02(IB02(ID_same_cs));

totalmzList=[];
for j=1:length(Pepcommon0102v1)
    [peptidenew, massdiffList, isHeavy]=modprocess({Pepcommon0102v1{j}});
    peptidenew=peptidenew{1};
    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(Pepcommon0102v1{j}));
    csvalue=Orbi_chargestate01v1(j);
    mzvalue=(weight+csvalue*1.0073)/csvalue;
    iso(j,:)=isotopepattern;
    totalmzList=[totalmzList; mzvalue mzvalue+1.00335/csvalue mzvalue+2.00547/csvalue mzvalue+3.00882/csvalue mzvalue+4.00849/csvalue mzvalue+5.01184/csvalue];
end
% finalmzList=[finalmzList finalmzList+1.00335./finalchargeList finalmzList+2.00547./finalchargeList finalmzList+3.00882./finalchargeList finalmzList+4.00849./finalchargeList finalmzList+5.01184./finalchargeList];

%%%%%%%%%% get XICs
totalmzList;
tolerance=20;
for i=1:6
    XICs_Orbi01{i}=getXICs(Orbi_data01l1,totalmzList(:,i),tolerance);
end
for i=1:6
    XICs_Orbi02{i}=getXICs(Orbi_data02l1,totalmzList(:,i),tolerance);
end
save XICs_Orbitrap_5569_f4 XICs_Orbi01 XICs_Orbi02

load XICs_Orbitrap_5569_f4
% for i=1:6
%     XICs_TOF03{i}=getXICs(Orbi_peakl03,totalmzList(:,i),tolerance);
% end
%%%%%%%%%%%

%%%%%%%%%%% interval detection
[IntervalList01,posi01]=Orbi_Verification_intervaldetection(Pepcommon0102v1,Orbi_ms2information01v1,XICs_Orbi01,Orbi_retentiont01l1);
[IntervalList02,posi02]=Orbi_Verification_intervaldetection(Pepcommon0102v1,Orbi_ms2information02v1,XICs_Orbi02,Orbi_retentiont02l1);
[V_p,IA,IB]=intersect(posi01,posi02);
        
Orbi_peplist01_final=Pepcommon0102v1(V_p);
Orbi_peplist02_final=Pepcommon0102v1(V_p);

num_c=0;  ID_c=[];
for i=1:length(Orbi_peplist01_final)-1

    String1=Orbi_peplist01_final{i};
    String2=Orbi_peplist02_final{i+1};    
    
    if length(String1)>=length(String2)
        Po_c=find(String1=='c');
        if ~isempty(Po_c)
            String1_cmp=String1(3:Po_c-1);
            String2_cmp=String2(3:end-2);
            JUD_cmp=strcmp(String1_cmp,String2_cmp);
            if JUD_cmp==1
                num_c=num_c+1;
                ID_c=[ID_c;i];
            end
        end
    else
        Po_c=find(String2=='c');
        if ~isempty(Po_c)
            String1_cmp=String1(3:end-2);
            String2_cmp=String2(3:Po_c-1);
            JUD_cmp=strcmp(String1_cmp,String2_cmp);
            if JUD_cmp==1
                num_c=num_c+1;
                ID_c=[ID_c;i+1];
            end
        end
    end
end
Orbi_peplist01_final(ID_c)=[];
Orbi_peplist02_final(ID_c)=[];
V_p_c=V_p;         V_p_c(ID_c)=[];
Orbi_protein01_final=Orbi_protein01v1(V_p_c);
Orbi_protein02_final=Orbi_protein02v1(V_p_c);
Orbi_chargestate01_final=Orbi_chargestate01v1(V_p_c);
Orbi_chargestate02_final=Orbi_chargestate02v1(V_p_c);
Orbi_ms2information01_final=Orbi_ms2information01v1(V_p_c);
Orbi_ms2information02_final=Orbi_ms2information02v1(V_p_c);
Orbi_IntervalList01_final=IntervalList01(V_p_c);
Orbi_IntervalList02_final=IntervalList02(V_p_c);
for i=1:6
    Orbi_XICs01{i}=XICs_Orbi01{i}(:,V_p_c);
    Orbi_XICs02{i}=XICs_Orbi02{i}(:,V_p_c);
end
Orbi_iso_final=iso(V_p_c,:);

totalmzList_Final=[];
for j=1:length(Orbi_peplist01_final)
    [peptidenew, massdiffList, isHeavy]=modprocess({Orbi_peplist01_final{j}});
    peptidenew=peptidenew{1};
    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(Orbi_peplist01_final{j}));
    csvalue=Orbi_chargestate01v1(j);
    mzvalue=(weight+csvalue*1.0073)/csvalue;
    totalmzList_Final=[totalmzList_Final; mzvalue mzvalue+1.00335/csvalue mzvalue+2.00547/csvalue mzvalue+3.00882/csvalue mzvalue+4.00849/csvalue mzvalue+5.01184/csvalue];
end

for i=1:6
    Orbi01_XICs01_final{i}=XICs_Orbi01{i}(:,V_p_c);
    Orbi02_XICs02_final{i}=XICs_Orbi02{i}(:,V_p_c);
end
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get ratio
for i=1:length(Orbi_peplist01_final)
    
    Interval_Scan_start=[];
    Interval_Scan_end=[];
    for nn=1:size(Orbi_IntervalList01_final{i}.Record,1)
        j=Orbi_IntervalList01_final{i}.Record(nn,1);
        k=Orbi_IntervalList01_final{i}.Record(nn,2);
        Interval_Scan_start=[Interval_Scan_start;Orbi_IntervalList01_final{i}.intervallist{j}(k,1)];
        Interval_Scan_end=[Interval_Scan_end;Orbi_IntervalList01_final{i}.intervallist{j}(k,2)];
    end    
    intervalsdata01=zeros(max(Interval_Scan_end)-min(Interval_Scan_start)+1,6);
    intervalsdata01v1=zeros(max(Interval_Scan_end)-min(Interval_Scan_start)+1,6);%%% cut through all the interval raw
    Intervaldata01_Scanstart=min(Interval_Scan_start);
    Intervaldata01_Scanend=max(Interval_Scan_end);

    Interval_Scan_start=[];
    Interval_Scan_end=[];
    for nn=1:size(Orbi_IntervalList02_final{i}.Record,1)
        j=Orbi_IntervalList02_final{i}.Record(nn,1);
        k=Orbi_IntervalList02_final{i}.Record(nn,2);
        Interval_Scan_start=[Interval_Scan_start;Orbi_IntervalList02_final{i}.intervallist{j}(k,1)];
        Interval_Scan_end=[Interval_Scan_end;Orbi_IntervalList02_final{i}.intervallist{j}(k,2)];
    end    
    intervalsdata02=zeros(max(Interval_Scan_end)-min(Interval_Scan_start)+1,6);
    intervalsdata02v1=zeros(max(Interval_Scan_end)-min(Interval_Scan_start)+1,6);%%% cut through all the interval raw
    Intervaldata02_Scanstart=min(Interval_Scan_start);
    Intervaldata02_Scanend=max(Interval_Scan_end);    
    
    for j=1:6
        
        intervalsdata01v1(:,j)=Orbi01_XICs01_final{j}(Intervaldata01_Scanstart:Intervaldata01_Scanend,i);
        intervalsdata02v1(:,j)=Orbi02_XICs02_final{j}(Intervaldata02_Scanstart:Intervaldata02_Scanend,i);
        Iso_posi=find(Orbi_IntervalList01_final{i}.Record(:,1)==j);
        if ~isempty(Iso_posi)
            k=Orbi_IntervalList01_final{i}.Record(Iso_posi,2);
            Scan_start01=Orbi_IntervalList01_final{i}.intervallist{j}(k,1);
            Scan_end01=Orbi_IntervalList01_final{i}.intervallist{j}(k,2);
            Interval_diff=Scan_end01-Scan_start01;
            intervalsdata01(Scan_start01-Intervaldata01_Scanstart+1:Scan_start01-Intervaldata01_Scanstart+1+Interval_diff,j)=Orbi01_XICs01_final{j}(Scan_start01:Scan_end01,i);
        end
        
        Iso_posi=find(Orbi_IntervalList02_final{i}.Record(:,1)==j);
        if ~isempty(Iso_posi)
            k=Orbi_IntervalList02_final{i}.Record(Iso_posi,2);
            Scan_start02=Orbi_IntervalList02_final{i}.intervallist{j}(k,1);
            Scan_end02=Orbi_IntervalList02_final{i}.intervallist{j}(k,2);
            Interval_diff=Scan_end02-Scan_start02;
            intervalsdata02(Scan_start02-Intervaldata02_Scanstart+1:Scan_start02-Intervaldata02_Scanstart+1+Interval_diff,j)=Orbi02_XICs02_final{j}(Scan_start02:Scan_end02,i);
        end       
        
    end
    
    %%%%%%%%%%%%%% interval data v1 is the all cut data 6 intervals based
    %%%%%%%%%%%%%% on maximum peak intensity;
        
    N_mean_filter=5;
    
    Intervaldata01_Scanstart;
    Intervaldata01_Scanend;
    
    intervalsdata01v1_afterfilter=zeros(size(intervalsdata01v1,1),size(intervalsdata01v1,2));
    for j=1:size(intervalsdata01v1,2)
        Weight=ones(N_mean_filter,1)./N_mean_filter;
        intervalsdata01v1_afterfilter(:,j)=filter2(Weight,intervalsdata01v1(:,j));    
    end
    
    [V01,P01]=min(abs(Orbi_retentiont01l1-Orbi_ms2information01_final(i)));
    for j=1:6
        MS2_Peak_Intensity01(j)=intervalsdata01v1(P01-Intervaldata01_Scanstart+1,j);
        MS2_Peak_Intensity01_afterfilter(j)=intervalsdata01v1_afterfilter(P01-Intervaldata01_Scanstart+1,j);
    end
    Spearman_Corr01=zeros(1,size(intervalsdata01v1,1));
    Spearman_Corr01_afterfilter=zeros(1,size(intervalsdata01v1,1));
    Top_isotope_num=6;
    for k=1:size(intervalsdata01v1,1)
        [V_ms2_Int01,P_ms2_Int01]=sort(MS2_Peak_Intensity01,'descend');
        a_toppeak=intervalsdata01v1(k,P_ms2_Int01(1:Top_isotope_num));
        Spearman_Corr01_toppeak(k)=corr(a_toppeak',MS2_Peak_Intensity01(P_ms2_Int01(1:Top_isotope_num))','type','spearman');
    
        a=intervalsdata01v1(k,:);
        Spearman_Corr01(k)=corr(a',MS2_Peak_Intensity01','type','spearman');    
        
        a_afterfilter=intervalsdata01v1_afterfilter(k,:);
        Spearman_Corr01_afterfilter(k)=corr(a_afterfilter',MS2_Peak_Intensity01_afterfilter','type','spearman'); 
    
    end
    
    Goodpoint_p01=Spearman_Corr01_afterfilter>=0.5;
%     Goodpoint_p01=Spearman_Corr01>=0.5;
%     Goodpoint_p01=Spearman_Corr01_toppeak>=0.5;
    
    Goodpoint_p01=[0,Goodpoint_p01,0];
    Rownumber01=P01-Intervaldata01_Scanstart+1+1;
    Rownumber01_start=Rownumber01;
    while Goodpoint_p01(Rownumber01_start)==1
        Rownumber01_start=Rownumber01_start-1;
    end
    Rownumber01_end=Rownumber01;
    while Goodpoint_p01(Rownumber01_end)==1
        Rownumber01_end=Rownumber01_end+1;
    end
    Goodinterval_start01=Rownumber01_start;
    Goodinterval_end01=Rownumber01_end-2;
    intervalsdata01v2=intervalsdata01v1(Goodinterval_start01:Goodinterval_end01,:);
   
    
    Intervaldata02_Scanstart;
    Intervaldata02_Scanend;
    intervalsdata02v1_afterfilter=zeros(size(intervalsdata02v1,1),size(intervalsdata02v1,2));
    for j=1:size(intervalsdata02v1,2)
        Weight=ones(N_mean_filter,1)./N_mean_filter;
        intervalsdata02v1_afterfilter(:,j)=filter2(Weight,intervalsdata02v1(:,j));    
    end
    
    [V02,P02]=min(abs(Orbi_retentiont02l1-Orbi_ms2information02_final(i)));
    for j=1:6
        MS2_Peak_Intensity02(j)=intervalsdata02v1(P02-Intervaldata02_Scanstart+1,j);
        MS2_Peak_Intensity02_afterfilter(j)=intervalsdata02v1_afterfilter(P02-Intervaldata02_Scanstart+1,j);
    end
    Spearman_Corr02=zeros(1,size(intervalsdata02v1,1));
    Spearman_Corr02_afterfilter=zeros(1,size(intervalsdata02v1,1));
    Top_isotope_num=6;
    for k=1:size(intervalsdata02v1,1)
        [V_ms2_Int02,P_ms2_Int02]=sort(MS2_Peak_Intensity02,'descend');
        a_toppeak=intervalsdata02v1(k,P_ms2_Int02(1:Top_isotope_num));
        Spearman_Corr02_toppeak(k)=corr(a_toppeak',MS2_Peak_Intensity02(P_ms2_Int02(1:Top_isotope_num))','type','spearman');
    
        a=intervalsdata02v1(k,:);
        Spearman_Corr02(k)=corr(a',MS2_Peak_Intensity02','type','spearman');    
        
        a_afterfilter=intervalsdata02v1_afterfilter(k,:);
        Spearman_Corr02_afterfilter(k)=corr(a_afterfilter',MS2_Peak_Intensity02_afterfilter','type','spearman'); 
    
    end
    
    Goodpoint_p02=Spearman_Corr02_afterfilter>=0.5;
%     Goodpoint_p02=Spearman_Corr02>=0.5;
%     Goodpoint_p02=Spearman_Corr021_toppeak>=0.5;

    Goodpoint_p02=[0,Goodpoint_p02,0];
    Rownumber02=P02-Intervaldata02_Scanstart+1+1;
    Rownumber02_start=Rownumber02;
    while Goodpoint_p02(Rownumber02_start)==1
        Rownumber02_start=Rownumber02_start-1;
    end
    Rownumber02_end=Rownumber02;
    while Goodpoint_p02(Rownumber02_end)==1
        Rownumber02_end=Rownumber02_end+1;
    end
    Goodinterval_start02=Rownumber02_start;
    Goodinterval_end02=Rownumber02_end-2;
    intervalsdata02v2=intervalsdata02v1(Goodinterval_start02:Goodinterval_end02,:);
    clear Spearman_Corr02
    
    %%%%%%%%%%%%%%
    %%%%% add zeros interval detection
    Orbi_data01_ms2Interval(i).intervalsdata=intervalsdata01;
    Orbi_data01_ms2Interval(i).iso=Orbi_iso_final(i,:);
    Orbi_data01_ms2Interval(i).Scanstart=Intervaldata01_Scanstart;
    Orbi_data01_ms2Interval(i).Scanend=Intervaldata01_Scanend;
    Orbi_data02_ms2Interval(i).intervalsdata=intervalsdata02;
    Orbi_data02_ms2Interval(i).iso=Orbi_iso_final(i,:);
    Orbi_data02_ms2Interval(i).Scanstart=Intervaldata02_Scanstart;
    Orbi_data02_ms2Interval(i).Scanend=Intervaldata02_Scanend;
    %%%%% based on maximum intensity interval detection
    Orbi_data01_ms2Intervalv1(i).intervalsdata=intervalsdata01v1;
    Orbi_data01_ms2Intervalv1(i).iso=Orbi_iso_final(i,:);
    Orbi_data01_ms2Intervalv1(i).Scanstart=Intervaldata01_Scanstart;
    Orbi_data01_ms2Intervalv1(i).Scanend=Intervaldata01_Scanend;
    Orbi_data02_ms2Intervalv1(i).intervalsdata=intervalsdata02v1;
    Orbi_data02_ms2Intervalv1(i).iso=Orbi_iso_final(i,:);
    Orbi_data02_ms2Intervalv1(i).Scanstart=Intervaldata02_Scanstart;
    Orbi_data02_ms2Intervalv1(i).Scanend=Intervaldata02_Scanend;
    %%%%% based on ms2 points spearman corr interval detection
    Orbi_data01_ms2Intervalv2(i).intervalsdata=intervalsdata01v2;
    Orbi_data01_ms2Intervalv2(i).iso=Orbi_iso_final(i,:);
    Orbi_data01_ms2Intervalv2(i).Scanstart=Intervaldata01_Scanstart+Goodinterval_start01-1;
    Orbi_data01_ms2Intervalv2(i).Scanend=Intervaldata01_Scanstart+Goodinterval_end01-1;
    Orbi_data02_ms2Intervalv2(i).intervalsdata=intervalsdata02v2;
    Orbi_data02_ms2Intervalv2(i).iso=Orbi_iso_final(i,:);
    Orbi_data02_ms2Intervalv2(i).Scanstart=Intervaldata02_Scanstart+Goodinterval_start02-1;
    Orbi_data02_ms2Intervalv2(i).Scanend=Intervaldata02_Scanstart+Goodinterval_end02-1;
    
    clear intervalsdata01 intervalsdata02
    clear intervalsdata01v1 intervalsdata02v1
    clear intervalsdata01v2 intervalsdata02v2
end

save newintervaldetect5569f4 Orbi_peplist01_final Orbi_peplist02_final Orbi_protein01_final Orbi_protein02_final Orbi_data01_ms2Intervalv2 Orbi_data02_ms2Intervalv2
