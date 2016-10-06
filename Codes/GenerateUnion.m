function [Totalinformationmatrix,Pepseq]=GenerateUnion(Total)
Num=1;
Totalinformationmatrix=[];
for i=1:length(Total.peptide)
    if Total.beta(i)==0
            Pepseq{Num}=Total.peptide{i};
            Num=Num+1;
            TF=strcmp(Total.peptide{i},Total.peptide);
            Posi=find(TF==1);
            inf_vector=zeros(1,24);
            if length(Posi)>1
                for k=1:length(Posi);
                    Index_data=Total.remarks(Posi(k));
                    inf_vector((Index_data-1)*8+1:Index_data*8)=[Total.assumed_charge(Posi(k)), Total.mass(Posi(k)), Total.timepoint(Posi(k)), Total.remarks(Posi(k)), Total.interval(Posi(k),:), Total.probability(Posi(k)), Total.Indexineachdata(Posi(k))];
                end        
                Totalinformationmatrix=[Totalinformationmatrix;inf_vector];
                Total.beta(Posi)=1;
            else
                Index_data=Total.remarks(Posi);
                inf_vector((Index_data-1)*8+1:Index_data*8)=[Total.assumed_charge(Posi), Total.mass(Posi), Total.timepoint(Posi), Total.remarks(Posi), Total.interval(Posi,:), Total.probability(Posi), Total.Indexineachdata(Posi)];
                Totalinformationmatrix=[Totalinformationmatrix;inf_vector];    
                Total.beta(Posi)=1;
            end    
    end
end
