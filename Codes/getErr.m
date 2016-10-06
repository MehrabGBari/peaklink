function err=getErr(intensity, p1,p2,p3, a1, rate, f)
a2=2*rate*f*(1-f)/(1+rate*(1-f)^2)*a1;
a3=rate*f*f/(1+rate*(1-f)^2)*a1;
mymean=a1*p1+a2*p2+a3*p3;
err=abs(intensity-mymean);
err(intensity>0)=err(intensity>0)./intensity(intensity>0);
err=mean(err);
