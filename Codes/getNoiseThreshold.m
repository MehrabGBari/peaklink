function noiseThreshold=getNoiseThreshold(smoothXICs,XICs,noiseThresholdLevel);
             noiseVector=smoothXICs-XICs;
            noiseMean=mean(noiseVector(noiseVector>0));
            noiseStd=std(noiseVector(noiseVector>0));
            noiseThreshold=noiseMean+noiseStd*noiseThresholdLevel;  
end