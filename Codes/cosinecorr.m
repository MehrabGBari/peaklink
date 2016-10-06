function corr=cosinecorr(vect1,vect2)
corr=-inf;
if sum(vect1)==0
    return;
end
if sum(vect2)==0
    return;
end
corr=sum(vect1.*vect2)/sqrt(sum(vect1.^2)*sum(vect2.^2));