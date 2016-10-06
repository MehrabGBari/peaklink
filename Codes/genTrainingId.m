function TrainingId = genTrainingId(Arry1, Arry2)

IDMwithMS2=[Arry1,Arry2];
IDVwithMS2=sum(IDMwithMS2,2);
TrainingId=find(IDVwithMS2==2);