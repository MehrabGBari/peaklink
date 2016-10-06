function [Cal_HLRatio_new,ID_95]=Outlayer_cutting(Cal_HLRatio,Std_num)


CaL_UP_bound=mean(log2(Cal_HLRatio(Cal_HLRatio~=0)))+Std_num*std(log2(Cal_HLRatio(Cal_HLRatio~=0)));
CaL_Down_bound=mean(log2(Cal_HLRatio(Cal_HLRatio~=0)))-Std_num*std(log2(Cal_HLRatio(Cal_HLRatio~=0)));

ID_95=find(log2(Cal_HLRatio)>=CaL_Down_bound & log2(Cal_HLRatio)<=CaL_UP_bound);
Cal_HLRatio_new=Cal_HLRatio(ID_95);
% mean(Cal_HLRatio_new)
% std(log2(Cal_HLRatio_new))
% length(Cal_HLRatio_new)
