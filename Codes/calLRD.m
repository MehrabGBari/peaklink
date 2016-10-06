function LRD = calLRD(uniqueProteinList, uniHlrMean)

[index1_r1, index2_r1] = uniquePeptideDetection(uniqueProteinList);
hlrG1tmp = uniHlrMean(index1_r1);
hlrG2tmp = uniHlrMean(index2_r1);
posNull = find(hlrG1tmp>0 & hlrG1tmp<1000 & hlrG2tmp > 0 & hlrG2tmp <1000);
hlrG1 = hlrG1tmp(posNull);
hlrG2 = hlrG2tmp(posNull);
LRD = std(log2(hlrG1) - log2(hlrG2));