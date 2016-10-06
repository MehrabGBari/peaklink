function [Percentage,STD_rand,MEAN_rand]=PERmuteTheRatio(Ratio)

num=1;
for P=0:0.04:1
    L=round(length(Ratio(:,1))*P);
    ID_picked=randsample(length(Ratio(:,1)),L,'false');
    V_01_pick=Ratio(ID_picked,1);
    V_02_pick=Ratio(ID_picked,2);
    Rand_V_01_pick=V_01_pick(randsample(L,L,'false'));
    Rand_V_02_pick=V_02_pick(randsample(L,L,'false'));
    Rand_Orbit_o18rate01=Ratio(:,1);
    Rand_Orbit_o18rate02=Ratio(:,2);
    Rand_Orbit_o18rate01(ID_picked)=Rand_V_01_pick;
    Rand_Orbit_o18rate02(ID_picked)=Rand_V_02_pick;
    STD_rand(num)=std(log2(Rand_Orbit_o18rate01)-log2(Rand_Orbit_o18rate02));
    MEAN_rand(num)=mean(log2(Rand_Orbit_o18rate01)-log2(Rand_Orbit_o18rate02));
    Percentage(num)=P;
    num=num+1;
end