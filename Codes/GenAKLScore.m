function [AKLScoreMatrixCorr,AKLScoreMatrixNCorr]=GenAKLScore(AKLmatrix,Para)

for i=1:size(AKLmatrix,1)
    for j=1:size(AKLmatrix,2)
        AKLScoreMatrixCorr(i,j)=normpdf(AKLmatrix(i,j),Para.Para_AKLmodel(1),Para.Para_AKLmodel(2));
        AKLScoreMatrixNCorr(i,j)=normpdf(AKLmatrix(i,j),Para.Para_AKLmodel(3),Para.Para_AKLmodel(4));        
    end
end
