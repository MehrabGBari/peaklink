function HLR = quantHLR(FinalResultSel)


HLR = zeros(length(FinalResultSel)  , 3);
for i=1:length(FinalResultSel)
    for j=1:3
        IntensitySum=sum(FinalResultSel{i}.ElutionProfile{j},1);
        if length(IntensitySum)==8
            HLR(i,j)=sum(sum(IntensitySum(5:6)))/sum(sum(IntensitySum(1:2)));
        else
            HLR(i,j)=0;
        end
    end
end