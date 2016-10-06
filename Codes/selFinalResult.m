% select final result


function [FinalResult_forTest, pos_forTest, posE] = selFinalResult(FinalResult, pos, flagA, flag)

posE = find(flagA(:, flag)==1);



FinalResult_forTest = FinalResult(posE);
pos_forTest = pos(posE);