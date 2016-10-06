function [SC10data,SC10result,Peptide_SC10_information]=Mascort_result_processing(Mascortfile_path)

[SC10data,SC10result]=readtext(Mascortfile_path,'\t');

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