function Iso12KL_th=Generate_iso12KL_th(iso)

num=1500;
RandID=randsample(1:size(iso,1),num,'false');
Average_iso(:,1)=iso(RandID,1);
Average_iso(:,2)=iso(RandID,2);
for i=1:size(Average_iso,1)
    Average_norm_iso(i,:)=Average_iso(i,:)/sum(Average_iso(i,:));
end

Range=100000;
num_generate=100000;
Generate_iso(:,1)=randsample(1:Range,num_generate)/Range;
Generate_iso(:,2)=randsample(1:Range,num_generate)/Range;
Norm_Generate_iso(:,1)=Generate_iso(:,1)./sum(Generate_iso,2);
Norm_Generate_iso(:,2)=Generate_iso(:,2)./sum(Generate_iso,2);
% Generate_iso(:,2)=1-Generate_iso(:,1);
nn=1;
for kk=1:5
    for i=1:size(Average_norm_iso,1)
        ID_rand=randint(1,1,[1 size(Norm_Generate_iso,1)]);
        log_KL_Value(nn)=KL_calculate(Norm_Generate_iso(ID_rand,:),Average_norm_iso(i,:));
        nn=nn+1;
    end
end

mean(log_KL_Value);
figure
hist(log_KL_Value,1000)

Sort_Log_KL=sort(log_KL_Value);
Sort_Log_KL(round(length(log_KL_Value)*0.05))

[Y,I]=hist(log_KL_Value,1000);
for i=1:1000
    PER(i)=sum(Y(1:i))/sum(Y);    
end
ID=find(PER>=0.05);
Iso12KL_th=I(ID(1));
figure
plot(PER)

