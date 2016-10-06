function [ATScoreMatrixCorr,ATScoreMatrixNCorr]=GenATScore(Tmatrix,Para)

for i=1:size(Tmatrix,1)
    for j=1:size(Tmatrix,2)
        ATScoreMatrixCorr(i,j)=normpdf(Tmatrix(i,j),Para.Para_ATmodel(1),Para.Para_ATmodel(2));
        ATScoreMatrixNCorr(i,j)=normpdf(Tmatrix(i,j),Para.Para_ATmodel(3),Para.Para_ATmodel(4));        
    end
end


