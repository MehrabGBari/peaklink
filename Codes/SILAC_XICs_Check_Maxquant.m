clc
clear all

[Maxquant_10data, Maxquant_10result]=readtext('C:\Users\Jian Cui\Desktop\Xiaolin Expt01212012\A+B10_Maxquant\combined\txt\peptides.txt','\t');
Maxquant_Razor_Protein=Maxquant_10data(2:end,34);
Maxquant_Razor_Protein_HLR=Maxquant_10data(2:end,43);
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
    Final_Maxquant_Razor_Protein_unique_HLR(i)=sum(Final_Maxquant_Razor_Protein_HLR(ID))/length(ID);
end
length(Final_Maxquant_Razor_Protein_unique_HLR)
mean(log2(Final_Maxquant_Razor_Protein_unique_HLR))
std(log2(Final_Maxquant_Razor_Protein_unique_HLR))

Maxquant_peptide=Maxquant_10data(2:end,8);
Maxquant_peptide_HLR=Maxquant_10data(2:end,43);
Maxquant_peptide_cs=Maxquant_10data(2:end,40);
Maxquant_peptide_mass=Maxquant_10data(2:end,32);

SALIC_totalmzList=[];
for i=1:length(Maxquant_peptide_mass)    
    
    Pep_seq=Maxquant_peptide{i};
    K_posi=find(Pep_seq=='K');
    R_posi=find(Pep_seq=='R');
%     Mod_inform=Peptide_SC10_information.pep_var_mod{i};
%     K_Mod_string='Label:13C(6) (K)';
%     R_Mod_string='Label:13C(6) (R)';
%     ID_K_Mod_str_match=strfind(Mod_inform,K_Mod_string);
%     ID_R_Mod_str_match=strfind(Mod_inform,R_Mod_string);    
    N_K=length(K_posi);
    N_R=length(R_posi);
    mzvalue=(Maxquant_peptide_mass{i}+Maxquant_peptide_cs{i}*1.0073)/Maxquant_peptide_cs{i};
    SALIC_totalmzList=[SALIC_totalmzList; mzvalue mzvalue+(13.0034-12)/csvalue mzvalue+1.0034*2/csvalue mzvalue+1.0034*3/csvalue mzvalue+(N_K+N_R)*6.020129/csvalue mzvalue+(13.0034-12)/csvalue+(N_K+N_R)*6.020129/csvalue mzvalue+1.0034*2/csvalue+(N_K+N_R)*6.020129/csvalue  mzvalue+1.0034*3/csvalue+(N_K+N_R)*6.020129/csvalue];
end

[AB10retentiontl1,AB10datal1,AB10peakl,AB10retentiont,AB10MZInt_l1l2]=readrawdata('C:\Users\Jian Cui\Desktop\Xiaolin Expt01212012\A+B10.mzXML');
save D:\Program\QTOF_replicate_identification\Xiaolin_Expt01212012\AB10information AB10retentiontl1 AB10datal1 AB10peakl AB10retentiont AB10MZInt_l1l2

load D:\Program\QTOF_replicate_identification\Xiaolin_Expt01212012\AB10information

tolerance=20;
for i=1:8
    SILAC_AB10_XICs_Tolerance{i}=getXIC_LC_new(AB10peakl,SALIC_totalmzList(:,i),tolerance);
end
save D:\Program\QTOF_replicate_identification\Xiaolin_Expt01212012\AB10XICs SILAC_AB10_XICs_Tolerance

size(SILAC_AB10_XICs_Tolerance{1})
colorarray=['r', 'k', 'g', 'b', 'm', 'y','c','r:'];
height=1000000;
for kkk=204:206%1:length(posi01)
    pep_seq_string=Maxquant_peptide{kkk};
    figure
    for i=1:8
        plot(AB10retentiontl1,SILAC_AB10_XICs_Tolerance{i}(:,kkk),colorarray(i))
        hold on
    end
    grid on
    
%     stem(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1)),height*2*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'r*')
%     stem(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)),height*2*ones(size(IntervalList01{kk}.intervallist_after7_combine,1),1),'k*')
%     X=(SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,1))+SC10retentiontl1(IntervalList01{kk}.intervallist_after7_combine(:,2)))/2;
%     for RRR=1:length(X)
%         text(X(RRR),height*3,num2str(IntervalList01{kk}.Labelling_efficiency_after7_combine(RRR)))
%     end
    title(pep_seq_string)

end







