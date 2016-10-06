clc
clear all

%%%%%%%%%%%%%%%%% Maxquant Result processing
%%%%%%%%%%%%%%%%% A+B(1:1)     BB(1:5)    CC(5:1)
%%%%%%%%%%%%%%%%% A+B[10(done),25(done),70(done)]; BB[25(done),35(done),50(done),70(done)]; CC[35(done),50(done),70(done)];

% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B10_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B25_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B70_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB25_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB35_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB50_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB70_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC35_data\combined\txt\peptides.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC50_data\combined\txt\peptides.txt';
Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC70_data\combined\txt\peptides.txt';

[Maxquant_data, Maxquant_result,...
    Final_Maxquant_Razor_Protein,...
    Final_Maxquant_Razor_Protein_HLR,...
    Final_Maxquant_Razor_Protein_unique_HLR]=Maxquant_result_processing(Maxquantfile_path);

Maxquant_peptide_sequence=Maxquant_data(2:end,9);
Maxquant_Razor_Protein_HLR=Maxquant_data(2:end,44);
Maxquant_HLRatio=[];
Klabeled_id=[];
for i=1:length(Maxquant_Razor_Protein_HLR)
    Pep_seq=Maxquant_peptide_sequence{i};
    K_posi=find(Pep_seq=='K');
    if ~isempty(K_posi)        
        Klabeled_id=[Klabeled_id; i];
    end
    if ~isnan(Maxquant_Razor_Protein_HLR{i})
        Maxquant_HLRatio=[Maxquant_HLRatio; Maxquant_Razor_Protein_HLR{i},i];        
    end    
end
% length(Maxquant_Razor_Protein_HLR)
% size(Maxquant_HLRatio)
% size(Klabeled_id)
% length(intersect(Klabeled_id,Maxquant_HLRatio(:,2)))

[Klabeled_id_inter,IA_KL,IB_Max_HLR]=intersect(Klabeled_id,Maxquant_HLRatio(:,2));

% size(Maxquant_HLRatio)
% length(unique(Maxquant_peptide_sequence))
% figure
% hist(log2(Maxquant_HLRatio(:,1)),[-20:0.1:10]);grid on
% Id_g=find(log2(Maxquant_HLRatio(:,1))>-1);
% length(Id_g)
% mean((Maxquant_HLRatio(:,1)))
% std(log2(Maxquant_HLRatio(:,1)))

% length(Final_Maxquant_Razor_Protein_unique_HLR)
% mean((Final_Maxquant_Razor_Protein_unique_HLR))
% std(log2(Final_Maxquant_Razor_Protein_unique_HLR))
% figure
% hist(log2(Final_Maxquant_Razor_Protein_unique_HLR),[-20:0.1:10]);grid on

%%%%%%%%%%%%% read msms txt file
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B10_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B25_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B70_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB25_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB35_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB50_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB70_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC35_data\combined\txt\msms.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC50_data\combined\txt\msms.txt';
Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC70_data\combined\txt\msms.txt';


[Maxquant_data_msms, Maxquant_result_msms]=readtext(Maxquantfile_path,'\t');
Maxquant_msmspep=Maxquant_data_msms(2:end,12);
Maxquant_msmspep_scannumber=Maxquant_data_msms(2:end,10);
Final_Maxquant_msmspep=unique(Maxquant_msmspep);
%%%% the Final_Maxquant_msmspep is listed in the file
%%%% Maxquant_data peptide

%%%%%%%%%%%% read Oxidation txt file
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B10_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B25_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B70_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB25_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB35_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB50_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB70_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC35_data\combined\txt\Oxidation (M)Sites.txt';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC50_data\combined\txt\Oxidation (M)Sites.txt';
Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC70_data\combined\txt\Oxidation (M)Sites.txt';

[Maxquant_data_Oxidation, Maxquant_result_Oxidation]=readtext(Maxquantfile_path,'\t');
Maxquant_data_Oxidation_pep_Mo=Maxquant_data_Oxidation(2:end,33);
for i=1:length(Maxquant_data_Oxidation_pep_Mo)
    ID=find(Maxquant_data_Oxidation_pep_Mo{i}>=65 & Maxquant_data_Oxidation_pep_Mo{i}<=90);
    Maxquant_data_Oxidation_pep{i}=Maxquant_data_Oxidation_pep_Mo{i}(ID);
end
Maxquant_data_Oxidation_mz=Maxquant_data_Oxidation(2:end,37);
%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%% Mascort Result processing
% Mascortfile_path='D:\Program\QTOF_replicate_identification\Xiaolin_Expt01122012_Mascort\sc10.txt';
% [SC10data,SC10result,Peptide_SC10_information]=Mascort_result_processing(Mascortfile_path);
% save D:\Program\QTOF_replicate_identification\haskin_new_salic_data\Peptide_SC10_information Peptide_SC10_information
% load D:\Program\QTOF_replicate_identification\haskin_new_salic_data\Peptide_SC10_information
% 
% Mascortfile_distriller_path='D:\Program\QTOF_replicate_identification\Xiaolin_Expt01122012_Mascort\A+B10.txt';
% [_dis_data,_dis_result,Protein_dis_information]=Mascort_distriller_result_processing(Mascortfile_distriller_path);
% %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Generate XICs based on peptide list (can be Mascort peptide list or Maxquant peptide list)
Maxquant_peptide=Maxquant_data(2:end,9);
Maxquant_peptide_Oxidation=Maxquant_data(2:end,8);
Maxquant_peptide_HLR=Maxquant_data(2:end,44);
Maxquant_peptide_cs=Maxquant_data(2:end,41);
Maxquant_peptide_mass=Maxquant_data(2:end,33);
M_Oxidation_id=[];
for i=1:length(Maxquant_peptide_Oxidation)
    Oxi_M_posi=Maxquant_peptide_Oxidation{i};
    if ~isempty(Oxi_M_posi)
         M_Oxidation_id=[M_Oxidation_id; i];
    end
end

for i=1:length(Maxquant_peptide)
%     [peptidenew, massdiffList, isHeavy]=modprocess({Maxquant_peptide{i}});
%     peptidenew=peptidenew{1};
%     [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,4,getmodificationformula(Maxquant_peptide{i}));
    [peptideformula,isotopepattern,weight]=aminocalculation(Maxquant_peptide{i},4);
    Maxquant_iso(i,:)=isotopepattern;
    Long_mass(i)=weight;
%     csvalue=AMT_database_5569_f4_inLongTOFpep(i).charge;
%     mzvalue=(weight+csvalue*1.0073)/csvalue;
%     TOF_totalmzList=[TOF_totalmzList; mzvalue mzvalue+1.00335/csvalue mzvalue+2.00547/csvalue mzvalue+3.00882/csvalue mzvalue+4.00849/csvalue mzvalue+5.01184/csvalue];
end

for i=1:length(Maxquant_peptide)
    BestMSMSID=Maxquant_data{i+1,6};
    Maxquant_peptide_ms2scannumber(i)=Maxquant_msmspep_scannumber{BestMSMSID+1};
end

SALIC_totalmzList=[];Oxida_Id=[];
K_labeled_ID=[];No_K_labeled_ID=[];
for i=1:length(Maxquant_peptide_mass)    
    
    Pep_seq=Maxquant_peptide{i};
    K_posi=find(Pep_seq=='K');
    R_posi=find(Pep_seq=='R');
    if ~isempty(K_posi)
        K_labeled_ID=[K_labeled_ID; i];
    else No_K_labeled_ID=[No_K_labeled_ID; i];
    end
    LID_M_oxidation=strcmp(Pep_seq,Maxquant_data_Oxidation_pep);
    ID_M_oxidation=find(LID_M_oxidation==1);
%     Mod_inform=Peptide_SC10_information.pep_var_mod{i};
%     K_Mod_string='Label:13C(6) (K)';
%     R_Mod_string='Label:13C(6) (R)';
%     ID_K_Mod_str_match=strfind(Mod_inform,K_Mod_string);
%     ID_R_Mod_str_match=strfind(Mod_inform,R_Mod_string);    
    N_K=length(K_posi);
%     N_R=length(R_posi);%%%%% R is labeled
    N_R=0;%%%%% R is not labeled
    csvalue=Maxquant_peptide_cs{i}(1);
    if ischar(csvalue)
        csvalue=str2num(csvalue);
    end
    if sum(LID_M_oxidation)~=0
        Oxida_Id=[Oxida_Id;i];
       mzvalue=Maxquant_data_Oxidation_mz{ID_M_oxidation};
    else
        mzvalue=(Maxquant_peptide_mass{i}+csvalue*1.0073)/csvalue;
    end
    SALIC_totalmzList=[SALIC_totalmzList; mzvalue mzvalue+(13.0034-12)/csvalue mzvalue+1.0034*2/csvalue mzvalue+1.0034*3/csvalue mzvalue+(N_K+N_R)*6.020129/csvalue mzvalue+(13.0034-12)/csvalue+(N_K+N_R)*6.020129/csvalue mzvalue+1.0034*2/csvalue+(N_K+N_R)*6.020129/csvalue  mzvalue+1.0034*3/csvalue+(N_K+N_R)*6.020129/csvalue];

end
%%%%%%%%% A+B(1:1)     BB(1:5)    CC(5:1)
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B10_data\A+B10.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B25_data\A+B25.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data02\A+B70_data\A+B70.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB25_data\BB25.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB35_data\BB35.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB50_data\BB50.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\BB70_data\BB70.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC35_data\CC35.mzXML';
% Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC50_data\CC50.mzXML';
Maxquantfile_path='C:\Users\Jian Cui\Desktop\Haskin_data03\CC70_data\CC70.mzXML';

% [retentiontl1,datal1,peakl,retentiont,MZInt_l1l2]=readrawdata(Maxquantfile_path);

% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B10_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B25_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B70_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB25_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB35_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB50_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB70_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC35_result';
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC50_result';
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC70_result';

% save([Result_Save_Path,'\TRY'], 'tolerance')
% save([Result_Save_Path,'\information'], 'retentiontl1', 'datal1', 'peakl', 'retentiont', 'MZInt_l1l2')
% load([Result_Save_Path,'\information']);

tolerance=20;
SALIC_totalmzList_K_labeled=SALIC_totalmzList(K_labeled_ID,:);
Maxquant_peptide_K_labeled=Maxquant_peptide(K_labeled_ID);
Maxquant_iso_K_labeled=Maxquant_iso(K_labeled_ID,:);
Maxquant_peptide_HLR_K_labeled=Maxquant_peptide_HLR(K_labeled_ID);
Maxquant_peptide_cs_K_labeled=Maxquant_peptide_cs(K_labeled_ID);
Maxquant_peptide_mass_K_labeled=Maxquant_peptide_mass(K_labeled_ID);
Maxquant_peptide_ms2scannumber_K_labeled=Maxquant_peptide_ms2scannumber(K_labeled_ID);

pep01=Maxquant_peptide_K_labeled;
save([Result_Save_Path,'\Pepinformation'], 'pep01');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PredefindeRatio=1.08;
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B10_result';
% [Peptideseq_AB10,Resultmatrix_AB10]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B25_result';
[Peptideseq_AB25,Resultmatrix_AB25]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B70_result';
[Peptideseq_AB70,Resultmatrix_AB70]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Pep_11=[Peptideseq_AB25,Peptideseq_AB70];
Information_M_11=[Resultmatrix_AB25;Resultmatrix_AB70];
PredefindeRatio=0.18;
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB25_result';
[Peptideseq_BB25,Resultmatrix_BB25]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB35_result';
[Peptideseq_BB35,Resultmatrix_BB35]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB50_result';
[Peptideseq_BB50,Resultmatrix_BB50]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB70_result';
[Peptideseq_BB70,Resultmatrix_BB70]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Pep_15=[Peptideseq_BB25,Peptideseq_BB35,Peptideseq_BB50,Peptideseq_BB70];
Information_M_15=[Resultmatrix_BB25;Resultmatrix_BB35;Resultmatrix_BB50;Resultmatrix_BB70];
PredefindeRatio=4.56;
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC35_result';
[Peptideseq_CC35,Resultmatrix_CC35]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC50_result';
[Peptideseq_CC50,Resultmatrix_CC50]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC70_result';
[Peptideseq_CC70,Resultmatrix_CC70]=Read_cal_Max_HLR_raw(Result_Save_Path,PredefindeRatio);
Pep_51=[Peptideseq_CC35,Peptideseq_CC50,Peptideseq_CC70];
Information_M_51=[Resultmatrix_CC35;Resultmatrix_CC50;Resultmatrix_CC70];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% % for i=1:8
% %     SILAC_XICs_Tolerance{i}=getXIC_LC_new(peakl,SALIC_totalmzList_K_labeled(:,i),tolerance);
% % end
% % save([Result_Save_Path,'\XICs_Klabeled'], 'SILAC_XICs_Tolerance');
% load([Result_Save_Path,'\XICs_Klabeled']);
% 
% % size(SILAC_XICs_Tolerance{1})
% % colorarray=['r', 'k', 'g', 'b', 'm', 'y','c','b:'];
% % height=1000000;
% % Maxquant_detectedHLR_ID=Maxquant_HLRatio(:,2);
% % Maxquant_undetectedHLR_ID=1:size(SILAC_XICs_Tolerance{1},2);
% % Maxquant_undetectedHLR_ID(Maxquant_detectedHLR_ID)=[];
% % for iii=101:105%1:length(posi01)
% %     kkk=Maxquant_undetectedHLR_ID(iii);
% % %     kkk=Maxquant_detectedHLR_ID(iii);
% %     pep_seq_string=Maxquant_peptide{kkk};
% %     figure
% %     for i=1:8
% %         plot(retentiontl1,SILAC_XICs_Tolerance{i}(:,kkk),colorarray(i))
% %         hold on
% %     end
% %     grid on
% %     
% % %     stem(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1)),height*2*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'r*')
% % %     stem(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)),height*2*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'k*')
% % %     X=(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1))+SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)))/2;
% % %     for RRR=1:length(X)
% % %         text(X(RRR),height*3,num2str(IntervalList01{kk}.Labelling_efficiency_after7_combine(RRR)))
% % %     end
% % 
% %     stem(retentiont(Maxquant_peptide_ms2scannumber(kkk),1),height*2,'r*')
% %     title(pep_seq_string)
% % 
% % end
% % 
% % Maxquant_iso(601:605,:);
% %%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%% New interval detection
% % for i=1:length(Peptide_SC10_information.pep_seq)    
% %     pep01{i}=[Peptide_SC10_information.pep_res_before{i},'.',Peptide_SC10_information.pep_seq{i},'.',Peptide_SC10_information.pep_res_after{i}];
% %     Pep_unique(i)=Peptide_SC10_information.pep_isunique{i};
% % end
% % ID_unique=find(Pep_unique==1);
% % length(unique(pep01(ID_unique)))
% pep01=Maxquant_peptide_K_labeled;
% posi01v1=[];
% posi01v2=[];
% posi01=[];
% IntervalList01=cell(1,length(pep01));
% for i=1:length(pep01)
%     SILAC_XICs01_sepc_pep=zeros(size(SILAC_XICs_Tolerance{1},1),8);
%     for j=1:8
%         SILAC_XICs01_sepc_pep(:,j)=SILAC_XICs_Tolerance{j}(:,i);        
%     end
%     [IntervalList,Jud_detect_good]=SILAC_Verification_intervaldetection_Spec_Pep(pep01{i},Maxquant_iso_K_labeled(i,:),SILAC_XICs01_sepc_pep);
%     IntervalList01{i}=IntervalList;
% %     if Jud_detect_good==1
% %         posi01v1=[posi01v1;i];
% %     end
%     if IntervalList01{i}.intervallist_after7_combine(1,1)~=0 && IntervalList01{i}.intervallist_after7_combine(1,2)~=0
%         posi01v2=[posi01v2;i];
%         retentiontl1v1=retentiontl1';
%         Pep_ms2time=retentiont(Maxquant_peptide_ms2scannumber_K_labeled(i),1);
%         T_m=retentiontl1v1(IntervalList01{i}.intervallist_after7_combine)-Pep_ms2time;
%         T_v=T_m(:,1).*T_m(:,2);
%         Id_ms2interval=find(T_v<=0);
%         if ~isempty(Id_ms2interval)
%             posi01=[posi01;i Id_ms2interval];
%         end            
%     end
% end
% length(pep01)
% length(posi01)
% length(posi01v2)
% 
% save ([Result_Save_Path,'\IntervalList01'], 'IntervalList01', 'posi01', 'posi01v2');
% 
% % colorarray=['r', 'k', 'g', 'b', 'm', 'y','c','k'];
% % height=1000000;
% % for kkk=151:155%1:length(posi01)
% % %     kkk=Maxquant_undetectedHLR_ID(iii);
% % %     kkk=Maxquant_detectedHLR_ID(iii);
% %     pep_seq_string=Maxquant_peptide_K_labeled{kkk};
% %     figure
% %     for i=1:8
% %         plot(retentiontl1,SILAC_XICs_Tolerance{i}(:,kkk),colorarray(i))
% %         hold on
% %     end
% %     grid on
% %     
% %     stem(retentiontl1(IntervalList01{kkk}.intervallist_after7_combine(:,1)),height*ones(size(IntervalList01{kkk}.intervallist_after7_combine,1),1),'ro')
% %     stem(retentiontl1(IntervalList01{kkk}.intervallist_after7_combine(:,2)),height*ones(size(IntervalList01{kkk}.intervallist_after7_combine,1),1),'ko')
% % %     X=(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1))+SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)))/2;
% % %     for RRR=1:length(X)
% % %         text(X(RRR),height*3,num2str(IntervalList01{kk}.Labelling_efficiency_after7_combine(RRR)))
% % %     end
% % 
% %     stem(retentiont(Maxquant_peptide_ms2scannumber_K_labeled(kkk),1),height*2,'r*')
% %     title(pep_seq_string)
% % 
% % end
% 
% %%%%%%%%%%%%%%%%%%% calculate H/L ratio based on our detection
% DIGTAL_P=2;
% for i=1:length(posi01)
%     kk=posi01(i,1);
%     id_int=posi01(i,2);    
%     Scan=IntervalList01{kk}.intervallist_after7_combine(id_int,:);
%     for j=1:8
%         Interval_Matrix(:,j)=SILAC_XICs_Tolerance{j}(Scan(1):Scan(2),kk);
%     end
% %     figure
% %     plot(Interval_Matrix)
%     Intensity_vector=sum(Interval_Matrix,1);
%     Interval_Matrix_w=Matrix_Weight_func(Interval_Matrix);
%     Intensity_vector_w=sum(Interval_Matrix_w,1);
%     Cal_HLRatio(i)=sum(Intensity_vector(5:5+DIGTAL_P-1))/sum(Intensity_vector(1:1+DIGTAL_P-1));
%     Cal_HLRatio_w(i)=sum(Intensity_vector_w(5:5+DIGTAL_P-1))/sum(Intensity_vector_w(1:1+DIGTAL_P-1));
%     clear Interval_Matrix Interval_Matrix_w
%     PepCal_Maxquant_peptide_HLR(i)=Maxquant_peptide_HLR{posi01(i,1)};
% end
% 
% save ([Result_Save_Path,'\Result_HLRatio_SumPo',num2str(DIGTAL_P)], 'Cal_HLRatio', 'Cal_HLRatio_w', 'PepCal_Maxquant_peptide_HLR');
% % [Inter_ID_kL,I_Klabeled_id,ID_posi01_id]=intersect(Klabeled_id,posi01(:,1));
% % length(ID_posi01_id)
% Cal_HLRatio_notinf=Cal_HLRatio(~isinf(Cal_HLRatio));
% [Cal_HLRatio_new,Cal_ID_95]=Outlayer_cutting(Cal_HLRatio_notinf);
% % figure
% % subplot(2,1,1)
% % hist(log2(Cal_HLRatio(Cal_HLRatio~=0)), -10:0.1:10);grid on;
% % subplot(2,1,2)
% % hist(log2(Cal_HLRatio_new), -10:0.1:10);grid on;
% 
% length(Cal_HLRatio)
% length(Cal_ID_95)
% mean(Cal_HLRatio_new)
% std(log2(Cal_HLRatio_new))
% 
% length(IB_Max_HLR)
% Maxquant_HLR_K_Labeled=Maxquant_HLRatio(IB_Max_HLR,1);
% [Maxquant_HLR_K_Labeled_new,Maxquant_HLR_K_Labeled_new_ID_95]=Outlayer_cutting(Maxquant_HLR_K_Labeled);
% % figure
% % subplot(2,1,1)
% % hist(log2(Maxquant_HLR_K_Labeled(Maxquant_HLR_K_Labeled~=0)), -10:0.1:10);grid on;
% % subplot(2,1,2)
% % hist(log2(Maxquant_HLR_K_Labeled_new), -10:0.1:10);grid on;
% 
% length(Maxquant_HLR_K_Labeled)
% length(Maxquant_HLR_K_Labeled_new_ID_95)
% mean(Maxquant_HLR_K_Labeled_new)
% std(log2(Maxquant_HLR_K_Labeled_new))
% 
% [I_COM,IA_COM,IB_COM]=intersect(posi01(Cal_ID_95,1),IA_KL(Maxquant_HLR_K_Labeled_new_ID_95));
% length(Cal_HLRatio_new(IA_COM))
% mean(Cal_HLRatio_new(IA_COM))
% std(log2(Cal_HLRatio_new(IA_COM)))
% length(Maxquant_HLR_K_Labeled_new(IB_COM))
% mean(Maxquant_HLR_K_Labeled_new(IB_COM))
% std(log2(Maxquant_HLR_K_Labeled_new(IB_COM)))
% 
% Cal_HLRatio_intersect_final=Cal_HLRatio_new(IA_COM);
% Maxquant_HLR_K_Labeled_intersect_final=Maxquant_HLR_K_Labeled_new(IB_COM);
% std(Cal_HLRatio_intersect_final)
% std(Maxquant_HLR_K_Labeled_intersect_final)
% save ([Result_Save_Path,'\Result_HLRatio_SumPo',num2str(DIGTAL_P),'_final'], 'Cal_HLRatio_intersect_final', 'Maxquant_HLR_K_Labeled_intersect_final');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B10_result';
% load([Result_Save_Path,'\IntervalList01']);
% load([Result_Save_Path,'\Result_HLRatio_SumPo',num2str(DIGTAL_P)]);
% 
% [Cal_HLRSumpo1_11_1,Maxquant_HLRSumpo1_11_1,...
%     Cal_HLRSumpo2_11_1,Maxquant_HLRSumpo2_11_1]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B25_result';
% [Cal_HLRSumpo1_11_2,Maxquant_HLRSumpo1_11_2,...
%     Cal_HLRSumpo2_11_2,Maxquant_HLRSumpo2_11_2]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data02\A+B70_result';
% [Cal_HLRSumpo1_11_3,Maxquant_HLRSumpo1_11_3,...
%     Cal_HLRSumpo2_11_3,Maxquant_HLRSumpo2_11_3]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB25_result';
% [Cal_HLRSumpo1_15_1,Maxquant_HLRSumpo1_15_1,...
%     Cal_HLRSumpo2_15_1,Maxquant_HLRSumpo2_15_1]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB35_result';
% [Cal_HLRSumpo1_15_2,Maxquant_HLRSumpo1_15_2,...
%     Cal_HLRSumpo2_15_2,Maxquant_HLRSumpo2_15_2]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB50_result';
% [Cal_HLRSumpo1_15_3,Maxquant_HLRSumpo1_15_3,...
%     Cal_HLRSumpo2_15_3,Maxquant_HLRSumpo2_15_3]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\BB70_result';
% [Cal_HLRSumpo1_15_4,Maxquant_HLRSumpo1_15_4,...
%     Cal_HLRSumpo2_15_4,Maxquant_HLRSumpo2_15_4]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC35_result';
% [Cal_HLRSumpo1_51_1,Maxquant_HLRSumpo1_51_1,...
%     Cal_HLRSumpo2_51_1,Maxquant_HLRSumpo2_51_1]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC50_result';
% [Cal_HLRSumpo1_51_2,Maxquant_HLRSumpo1_51_2,...
%     Cal_HLRSumpo2_51_2,Maxquant_HLRSumpo2_51_2]=Read_cal_Max_HLR(Result_Save_Path);
% Result_Save_Path='D:\Program\QTOF_replicate_identification\Haskins_data03\CC70_result';
% [Cal_HLRSumpo1_51_3,Maxquant_HLRSumpo1_51_3,...
%     Cal_HLRSumpo2_51_3,Maxquant_HLRSumpo2_51_3]=Read_cal_Max_HLR(Result_Save_Path);
% 
% 
% Pre_Value=4.56;
% Pred_v='5:1';
% Sum_Iso_num=1;
% 
% RA=[Pred_v(1),Pred_v(3)];
% eval(['AA=Cal_HLRSumpo', num2str(Sum_Iso_num), '_', num2str(RA), '_1;']);
% eval(['BB=Cal_HLRSumpo', num2str(Sum_Iso_num), '_', num2str(RA), '_2;']);
% eval(['CC=Cal_HLRSumpo', num2str(Sum_Iso_num), '_', num2str(RA), '_3;']);
% eval(['DD1=Maxquant_HLRSumpo', num2str(Sum_Iso_num), '_', num2str(RA), '_1;']);
% eval(['EE1=Maxquant_HLRSumpo', num2str(Sum_Iso_num), '_', num2str(RA), '_2;']);
% eval(['FF1=Maxquant_HLRSumpo', num2str(Sum_Iso_num), '_', num2str(RA), '_3;']);
% 
% DD=DD1';EE=EE1';FF=FF1';
% 
% for i=1:length(AA)
% %     Origin_name_D1S1{i}=['Data 1 Sum ',num2str(Sum_Iso_num),' iso'];
%     Origin_name_D1S1{i}='Proposed 1';
% end
% for i=1:length(BB)
% %     Origin_name_D2S1{i}=['Data 2 Sum ',num2str(Sum_Iso_num),' iso'];
%     Origin_name_D2S1{i}='Proposed 2';
% end
% for i=1:length(CC)
% %     Origin_name_D3S1{i}=['Data 3 Sum ',num2str(Sum_Iso_num),' iso'];
%     Origin_name_D3S1{i}='Proposed 3';
% end
% for i=1:length(DD)
% %     Origin_name_D1M{i}='Data 1 Maxquant';
%     Origin_name_D1M{i}='Maxquant 1';
% end
% for i=1:length(EE)
% %     Origin_name_D2M{i}='Data 2 Maxquant';
%     Origin_name_D2M{i}='Maxquant 2';
% end
% for i=1:length(FF)
% %     Origin_name_D3M{i}='Data 3 Maxquant';
%     Origin_name_D3M{i}='Maxquant 3';
% end
% 
% CMHLRXXXX=[AA,BB,CC,DD,EE,FF];
% AX=[Origin_name_D1S1,Origin_name_D2S1,Origin_name_D3S1,Origin_name_D1M,Origin_name_D2M,Origin_name_D3M];
% 
% HLR1=[AA,DD];
% AX1=[Origin_name_D1S1,Origin_name_D1M];
% HLR2=[BB,EE];
% AX2=[Origin_name_D2S1,Origin_name_D2M];
% HLR3=[CC,FF];
% AX3=[Origin_name_D3S1,Origin_name_D3M];
% 
% % figure
% % subplot(1,3,1)
% % boxplot(HLR1,AX1,'notch','on')
% % hold on
% % plot([0,1,2,3],[Pre_Value,Pre_Value,Pre_Value,Pre_Value],'k')
% % ylabel('H/L Ratio')
% % xlabel('Data 1')
% % subplot(1,3,2)
% % boxplot(HLR2,AX2,'notch','on')
% % hold on
% % plot([0,1,2,3],[Pre_Value,Pre_Value,Pre_Value,Pre_Value],'k')
% % ylabel('H/L Ratio')
% % xlabel('Data 2')
% % title(['Group (', Pred_v,') (Sum ',num2str(Sum_Iso_num),' isotope)'])
% % subplot(1,3,3)
% % boxplot(HLR3,AX3,'notch','on')
% % hold on
% % plot([0,1,2,3],[Pre_Value,Pre_Value,Pre_Value,Pre_Value],'k')
% % ylabel('H/L Ratio')
% % xlabel('Data 3')
% 
% figure
% boxplot(CMHLRXXXX,AX,'notch','on','widths',0.3)
% hold on
% plot([0,1,2,3,4,5,6,7],[Pre_Value,Pre_Value,Pre_Value,Pre_Value,Pre_Value,Pre_Value,Pre_Value,Pre_Value],'k')
% title(['Pre-calculate H/L ratio = ', num2str(Pre_Value)])
% ylabel('H/L Ratio')
% legend('Pre-calculate H/L ratio')
% 
% clear Origin_name_D1S1 Origin_name_D1M Origin_name_D2S1 Origin_name_D2M Origin_name_D3S1 Origin_name_D3M
% 
% 
% %%%%%%%%% plot sum result comparison
% Datax=[1 2 3 4];
% Data1y=[1.08 1.08 1.08 1.08];
% Data1y1=[1.0746 1.0097 0.9846 0.9800];
% Data1y2=[1.0794 1.0139 0.9877 0.9754];
% Data1y3=[1.1644 1.0801 1.0359 1.0313];
% Data1yy1=[0.4387 0.4313 0.4386 0.4531];
% Data1yy2=[0.4434 0.4231 0.4253 0.4276];
% Data1yy3=[0.4831 0.4665 0.4448 0.4690];
% 
% figure
% subplot(3,1,1)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y1, Datax, Data1yy1,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number') 
% title('Pre-calculate H/L Ratio=1.08');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% subplot(3,1,2)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y2, Datax, Data1yy2,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number');
% % title('Data 1:1 1.08');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% subplot(3,1,3)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y3, Datax, Data1yy3,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number') 
% % title('Data 1:1 1.08');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% 
% 
% Datax=[1 2 3 4];
% Data1y=[0.18 0.18 0.18 0.18];
% Data1y1=[0.2108 0.1873 0.1834 0.1852];
% Data1y2=[0.2141 0.1940 0.1906 0.1926];
% Data1y3=[0.2155 0.1996 0.1948 0.1933];
% Data1yy1=[0.4785 0.3895 0.4009 0.4263];
% Data1yy2=[0.4929 0.4492 0.4485 0.4886];
% Data1yy3=[0.5549 0.5194 0.5068 0.5011];
% 
% figure
% subplot(3,1,1)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y1, Datax, Data1yy1,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number') 
% title('Pre-calculate H/L Ratio=0.18');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% subplot(3,1,2)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y2, Datax, Data1yy2,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number');
% % title('Data 1:1 1.08');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% subplot(3,1,3)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y3, Datax, Data1yy3,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number') 
% % title('Data 1:1 1.08');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% 
% Datax=[1 2 3 4];
% Data1y=[4.56 4.56 4.56 4.56];
% Data1y1=[4.4069 4.1012 3.8734 3.7218];
% Data1y2=[4.2652 3.9575 3.7995 3.6904];
% Data1y3=[4.3439 4.0722 3.9174 3.7956];
% Data1yy1=[0.4314 0.3797 0.3780 0.4091];
% Data1yy2=[0.4855 0.4042 0.4301 0.4254];
% Data1yy3=[0.3358 0.2903 0.2818 0.2852];
% 
% figure
% subplot(3,1,1)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y1, Datax, Data1yy1,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% set(H1,'LineStyle','--','Marker','o');
% set(H2,'LineStyle',':','Marker','*');
% xlabel('Sum number')
% title('Pre-calculate H/L Ratio=4.56');
% subplot(3,1,2)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y2, Datax, Data1yy2,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number');
% % title('Data 1:1 1.08');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% subplot(3,1,3)
% plot(Datax, Data1y,'k-');legend('Pre-calculate ratio');hold on
% [AX,H1,H2] = plotyy(Datax, Data1y3, Datax, Data1yy3,'plot');
% set(get(AX(1),'Ylabel'),'String','Mean value') 
% set(get(AX(2),'Ylabel'),'String','Std value') 
% xlabel('Sum number') 
% % title('Data 1:1 1.08');
% set(H1,'LineStyle','--','Marker','o');set(H2,'LineStyle',':','Marker','*')
% 
% 
% 
% mean(log2(Cal_HLRatio(Cal_HLRatio~=0)))
% std(log2(Cal_HLRatio(Cal_HLRatio~=0)))
% mean(Cal_HLRatio(Cal_HLRatio<=20))
% std(log2(Cal_HLRatio(Cal_HLRatio<=20)))
% % mean((Cal_HLRatio_w))
% % std(log2(Cal_HLRatio_w))
% sum(Maxquant_HLRatio(:,1)<=10)
% length(Maxquant_HLRatio(:,1))
% mean((Maxquant_HLRatio(:,1)))
% std(log2(Maxquant_HLRatio(:,1)))
% 
% mean((Maxquant_HLR_K_Labeled))
% std(log2(Maxquant_HLR_K_Labeled))
% figure;hist(log2(Maxquant_HLR_K_Labeled),[-10:0.01:10])
% 
% [Inter_ID,I_Klabeled,ID_posi01]=intersect(Klabeled_id_inter,posi01(:,1));
% for i=1:length(Inter_ID)
%     KL_Maxquant_peptide_HLR(i)=Maxquant_Razor_Protein_HLR{Inter_ID(i),1};
% end
% Cal_HLRatio;
% figure
% subplot(2,1,1)
% hist(log2(Cal_HLRatio_w(ID_posi01)),[-20:0.1:5]);grid on
% subplot(2,1,2)
% hist(log2(KL_Maxquant_peptide_HLR),[-20:0.1:5]);grid on
% mean(log2(Cal_HLRatio_w(ID_posi01)))
% std(log2(Cal_HLRatio_w(ID_posi01)))
% mean((Cal_HLRatio(ID_posi01)))
% std(log2(Cal_HLRatio(ID_posi01)))
% mean((KL_Maxquant_peptide_HLR))
% std(log2(KL_Maxquant_peptide_HLR))
% 
% length(Cal_HLRatio)
% length(Maxquant_HLRatio(:,1))
% 
% L_ID_cal_Bzero=Cal_HLRatio>0;
% ID_cal_Bzero=find(L_ID_cal_Bzero==1);
% PepCal_Maxquant_peptide_HLR;
% L_ID_Maxquant_notnan=~isnan(PepCal_Maxquant_peptide_HLR);
% ID_Maxquant_notnan=find(L_ID_Maxquant_notnan==1);
% ID_Intersect=intersect(ID_Maxquant_notnan,ID_cal_Bzero);
% figure
% subplot(2,1,1)
% hist(log2(Cal_HLRatio(ID_Intersect)),[-20:0.1:5]);grid on
% subplot(2,1,2)
% hist(log2(PepCal_Maxquant_peptide_HLR(ID_Intersect)),[-20:0.1:5]);grid on
% ID_cal=log2(Cal_HLRatio(ID_Intersect))>=-1;
% ID_Maxquant=(log2(PepCal_Maxquant_peptide_HLR(ID_Intersect))>=-1);
% sum(log2(Cal_HLRatio(ID_Intersect))>=-1)
% sum(log2(PepCal_Maxquant_peptide_HLR(ID_Intersect))>=-1)
% mean((Cal_HLRatio(ID_Intersect(ID_cal))))
% std((Cal_HLRatio(ID_Intersect(ID_cal))))
% mean((PepCal_Maxquant_peptide_HLR(ID_Intersect(ID_Maxquant))))
% std((PepCal_Maxquant_peptide_HLR(ID_Intersect(ID_Maxquant))))
% 
% 
% 
% colorarray=['r', 'k', 'g', 'b', 'm', 'y','c','r:'];
% height=1000000;
% for kkk=1:length(posi01)
%     kk=posi01(kkk,1);
%     pep_seq_string=Maxquant_peptide{kk};
%     figure
%     for i=1:8
%         plot(retentiontl1,SILAC_XICs_Tolerance{i}(:,kk),colorarray(i))
%         hold on
%     end
%     grid on
%     
%     stem(retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1)),height*1*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'r*')
%     stem(retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)),height*1*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'k*')
%     X=(retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1))+retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)))/2;
%     for RRR=1:length(X)
%         text(X(RRR),height*3,num2str(IntervalList01{kk}.Labelling_efficiency_after7_combine(RRR)))
%     end
%     title(pep_seq_string)
%     
%     stem(retentiont(Maxquant_peptide_ms2scannumber(kk),1),height*2,'ro')
% 
% end
% %%%%%%%%%%%%%%%%%%%%










