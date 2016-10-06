function Information_Matrix_Total=CombineOriginalInformationMatrix(Information_Matrix_A,Information_Matrix_B)

[a,b]=size(Information_Matrix_A);
[c,d]=size(Information_Matrix_B);

pep_A=Information_Matrix_A(:,1);
pep_B=Information_Matrix_B(:,1);
ID_Same=[];
ID_notSame=[];
for i=1:length(pep_B)
    TF=strcmp(pep_B{i},pep_A);
    if sum(TF)>=1
        ID_inA=find(TF==1);
        ID_Same=[ID_Same; i ID_inA];
    else
        ID_notSame=[ID_notSame; i];        
    end    
end
Upright_matrix=num2cell(zeros(a,d));
Upright_matrix(ID_Same(:,2),:)=Information_Matrix_B(ID_Same(:,1),:);
downright_matrix=Information_Matrix_B(ID_notSame,:);
[e,f]=size(downright_matrix);
downleft_matrix=num2cell(zeros(e,b));
downleft_matrix(:,1)=downright_matrix(:,1);
Information_Matrix_Total=[Information_Matrix_A Upright_matrix; downleft_matrix downright_matrix];
