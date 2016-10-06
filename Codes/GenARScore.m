function [ARScoreMatrixCorr,ARScoreMatrixNCorr]=GenARScore(ARmatrix,Para)

for i=1:size(ARmatrix,1)
    for j=1:size(ARmatrix,2)
        ARScoreMatrixCorr(i,j)=gampdf(ARmatrix(i,j),Para.Para_ARmodel(1),Para.Para_ARmodel(2));
        ARScoreMatrixNCorr(i,j)=gampdf(ARmatrix(i,j),Para.Para_ARmodel(3),Para.Para_ARmodel(4));        
    end
end