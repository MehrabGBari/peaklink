function [LRD12, LRD23] = evaluateResult(HLR)
HLR_Vector=HLR(:,1);
HLR_Vectorv1=HLR_Vector(HLR_Vector>=0 & HLR_Vector<=1000);
figure
hist(HLR_Vectorv1,[0:0.1:20])

pos = find(HLR(:, 1) >0 & HLR(:, 1) <1000 & HLR(:, 2) >0 & HLR(:, 2) <1000  & HLR(:, 3) >0 & HLR(:, 3) <1000);
HLR_filtered = HLR(pos, :);
LRD12 = std(log2(HLR_filtered(:,1)) - log2(HLR_filtered(:,3)));
LRD23 = std(log2(HLR_filtered(:,2)) - log2(HLR_filtered(:,3)));
length(pos)





