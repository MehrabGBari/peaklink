function [Peplist_G1_com,ID_G1,ID_G1_Total]=Find_common_Pep(Peplist_G1,Final_peptide_list)
% Peplist_G1_com
% ID_G1
% ID_G1_Total
num=1;
for i=1:length(Peplist_G1)
    ID_com=strcmp(Peplist_G1{i}(3:end-2),Final_peptide_list);
    if sum(ID_com)~=0
        Peplist_G1_com{num}=Peplist_G1{i}(3:end-2);
        ID_G1(num)=i;
        ID_G1_Total(num)=find(ID_com==1);
        num=num+1;
    end    
end