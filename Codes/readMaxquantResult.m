function [Information_Matrix_K_labeled, Information_Matrix]=readMaxquantResult(path, PEPTIDE, MSMS, OXIDATION, LABEL_TECH)


%%%%%%%%%%%%%%%%%%%%%%%%% read in maxquant txt result
Maxquantfile_path=[path '\peptides.txt'];
[Maxquant_data, Maxquant_result,...
    Final_Maxquant_Razor_Protein,...
    Final_Maxquant_Razor_Protein_HLR,...
    Final_Maxquant_Razor_Protein_unique_HLR]=Maxquant_result_processing(Maxquantfile_path,PEPTIDE);

Maxquantfile_path=[path '\msms.txt'];
[Maxquant_data_msms, Maxquant_result_msms]=readtext(Maxquantfile_path,'\t');
Maxquant_msmspep=Maxquant_data_msms(2:end,MSMS.SEQ);%%%%%
Maxquant_msmspep_scannumber=Maxquant_data_msms(2:end,MSMS.SCANNUM);%%%%%
Final_Maxquant_msmspep=unique(Maxquant_msmspep);

Maxquantfile_path=[path '\Oxidation (M)Sites.txt'];

[Maxquant_data_Oxidation, Maxquant_result_Oxidation]=readtext(Maxquantfile_path,'\t');
Maxquant_data_Oxidation_pep_Mo=Maxquant_data_Oxidation(2:end,OXIDATION.OXIDATIONPROB);%%%%%
for i=1:length(Maxquant_data_Oxidation_pep_Mo)
    ID=find(Maxquant_data_Oxidation_pep_Mo{i}>=65 & Maxquant_data_Oxidation_pep_Mo{i}<=90);
    Maxquant_data_Oxidation_pep{i}=Maxquant_data_Oxidation_pep_Mo{i}(ID);
end
Maxquant_data_Oxidation_mz=Maxquant_data_Oxidation(2:end,OXIDATION.MZ);%%%%%

Maxquant_peptide=Maxquant_data(2:end,PEPTIDE.SEQ);%%%%%
Maxquant_peptide_Oxidation=Maxquant_data(2:end,PEPTIDE.OXIDATION);%%%%%

%Start Modification
Maxquant_peptide_HLR=ones(length(Maxquant_peptide),1);
Maxquant_peptide_HLR=mat2cell(Maxquant_peptide_HLR,ones(size(Maxquant_peptide_HLR,1),1),1);

%Maxquant_peptide_HLR=1;%Maxquant_data(2:end,PEPTIDE.HLR);%%%%%
% End Modification

Maxquant_peptide_cs=Maxquant_data(2:end,PEPTIDE.CS);%%%%%
Maxquant_peptide_mass=Maxquant_data(2:end,PEPTIDE.MASS);%%%%%
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
    BestMSMSID=Maxquant_data{i+1,PEPTIDE.BESTMSMSIDS};%%%%%
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
    
    if LABEL_TECH == 1
        N_R = 0;
    else if LABEL_TECH == 2
            N_R = length(R_posi);
        end
    end
%     if N_K~=0 && N_R ~= 0
%         i
%     end
%     N_R=length(R_posi);%%%%% R is labeled
    tmpA = Maxquant_peptide_cs{i};
    
    switch class(tmpA)      
        case 'double'  
            csvalue=Maxquant_peptide_cs{i};
        case 'char'
            pos = strfind(Maxquant_peptide_cs{i}, ',');
            if ~isempty(pos)
                csvalue = str2double(Maxquant_peptide_cs{i}(pos(1)-1));
            else
                csvalue=str2double(Maxquant_peptide_cs{i});
            end
    end




    if ischar(csvalue)
        csvalue=str2double(csvalue);
    end
    if sum(LID_M_oxidation)~=0
        Oxida_Id=[Oxida_Id;i];
       mzvalue=Maxquant_data_Oxidation_mz{ID_M_oxidation};
    else
        mzvalue=(Maxquant_peptide_mass{i}+csvalue*1.0073)/csvalue;
    end
    
    if LABEL_TECH == 1
        SALIC_totalmzList=[SALIC_totalmzList; ...
            mzvalue...
            mzvalue+(13.0034-12)/csvalue...
            mzvalue+1.0034*2/csvalue...
            mzvalue+1.0034*3/csvalue...
            mzvalue+(N_K+N_R)*6.020129/csvalue...
            mzvalue+(13.0034-12)/csvalue+(N_K+N_R)*6.020129/csvalue...
            mzvalue+1.0034*2/csvalue+(N_K+N_R)*6.020129/csvalue...
            mzvalue+1.0034*3/csvalue+(N_K+N_R)*6.020129/csvalue];
    else if LABEL_TECH == 2
              SALIC_totalmzList=[SALIC_totalmzList; ...
                  mzvalue...
                  mzvalue+(13.0034-12)/csvalue...
                  mzvalue+1.0034*2/csvalue...
                  mzvalue+1.0034*3/csvalue...
                  mzvalue+(N_K*8.0142+N_R*10.00827)/csvalue...
                  mzvalue+(13.0034-12)/csvalue+(N_K*8.0142+N_R*10.00827)/csvalue...
                  mzvalue+1.0034*2/csvalue+(N_K*8.0142+N_R*10.00827)/csvalue...
                  mzvalue+1.0034*3/csvalue+(N_K*8.0142+N_R*10.00827)/csvalue];
        end
    end
    csvalueFin(i) = csvalue;
end

%%%%%%%%%%%%%%% K-labeled peptide information

Information_Matrix_K_labeled(:,1)=Maxquant_peptide(K_labeled_ID);
Information_Matrix_K_labeled(:,2)=num2cell(csvalueFin(K_labeled_ID));
Information_Matrix_K_labeled(:,3)=Maxquant_peptide_mass(K_labeled_ID);
Information_Matrix_K_labeled(:,4)=num2cell(Maxquant_peptide_ms2scannumber(K_labeled_ID)');
Information_Matrix_K_labeled(:,5)=Maxquant_peptide_HLR(K_labeled_ID);
Information_Matrix_K_labeled(:,6)=Maxquant_data(K_labeled_ID+1,PEPTIDE.PROTEIN);
Information_Matrix_K_labeled(:,7)=num2cell(SALIC_totalmzList(K_labeled_ID,:),2);
Information_Matrix_K_labeled(:,8)=num2cell(Maxquant_iso(K_labeled_ID,:),2);

%%%%%%%%%%%%%%% Total peptide information
%%%%%%%%%%% column name: peptide_seq cs mass  ms2scannumber HLR protein totalmzlist iso

Information_Matrix(:,1)=Maxquant_peptide;
Information_Matrix(:,2)=num2cell(csvalueFin);
Information_Matrix(:,3)=Maxquant_peptide_mass;
Information_Matrix(:,4)=num2cell(Maxquant_peptide_ms2scannumber');
Information_Matrix(:,5)=Maxquant_peptide_HLR;
Information_Matrix(:,6)=Maxquant_data(2:end,PEPTIDE.PROTEIN);
Information_Matrix(:,7)=num2cell(SALIC_totalmzList,2);
Information_Matrix(:,8)=num2cell(Maxquant_iso,2);
