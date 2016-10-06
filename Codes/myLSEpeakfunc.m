function s=myLSEpeakfunc(a)
global p1;
global p2;
global p3;
%global sigma1;
%global sigma2;
%global sigma3;
%global ns;
global yintensity;
%a1=a(1);a2=a(2);a3=a(3);
a1=a(1);rate=a(2);f=a(3);
a2=2*rate*f*(1-f)/(1+rate*(1-f)^2)*a1;
a3=rate*f*f/(1+rate*(1-f)^2)*a1;
%myvar=a1*sigma1+a2*sigma2+a3*sigma3+ns;
mymean=a1*p1+a2*p2+a3*p3;
[temp1, idx1]=max([yintensity(1),yintensity(2)]);
[temp1, idx2]=max([yintensity(3),yintensity(4)]);
[temp1, idx3]=max([yintensity(5),yintensity(6)]);
peakidx=[idx1,2+idx2,4+idx3];
tempyintensity=yintensity(peakidx);
tempmymean=mymean(peakidx);
s=(tempyintensity-tempmymean)*(tempyintensity'-tempmymean');
end