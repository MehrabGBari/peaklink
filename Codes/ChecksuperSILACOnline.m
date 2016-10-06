clc
clear all

load SSMvsSuperSLIACinforamtion

ReplicateInformation_matrix=data(:,7:11);

Replicate01_Result=ReplicateInformation_matrix(:,1);
Replicate01_notnan=isnan(Replicate01_Result)~=1;
Replicate02_Result=ReplicateInformation_matrix(:,2);
Replicate02_notnan=isnan(Replicate02_Result)~=1;
Replicate03_Result=ReplicateInformation_matrix(:,3);
Replicate03_notnan=isnan(Replicate03_Result)~=1;
Replicate04_Result=ReplicateInformation_matrix(:,4);
Replicate04_notnan=isnan(Replicate04_Result)~=1;
Replicate05_Result=ReplicateInformation_matrix(:,5);
Replicate05_notnan=isnan(Replicate05_Result)~=1;

sum(Replicate01_notnan)
sum(Replicate02_notnan)
sum(Replicate03_notnan)
sum(Replicate04_notnan)
sum(Replicate05_notnan)
    
Total=Replicate01_notnan+Replicate02_notnan+Replicate03_notnan+Replicate04_notnan+Replicate05_notnan;
Zerosposi=Total==0; 
sum(Zerosposi)
    
sum(Replicate01_notnan)/5409
    
    
    
    
    
    
    