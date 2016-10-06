function [Maxquant_SC10data, Maxquant_SC10result,...
    Final_Maxquant_Razor_Protein,...
    Final_Maxquant_Razor_Protein_HLR,...
    Final_Maxquant_Razor_Protein_unique_HLR]=Maxquant_result_processing(Maxquantfile_path, PEPTIDE)


[Maxquant_SC10data, Maxquant_SC10result]=readtext(Maxquantfile_path,'\t');
Maxquant_Razor_Protein=Maxquant_SC10data(2:end,PEPTIDE.PROTEIN);
Maxquant_Razor_Protein_HLR=Maxquant_SC10data(2:end,PEPTIDE.SCORE);
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

% length(Final_Maxquant_Razor_Protein_unique_HLR)
% mean(log2(Final_Maxquant_Razor_Protein_unique_HLR))
% std(log2(Final_Maxquant_Razor_Protein_unique_HLR))