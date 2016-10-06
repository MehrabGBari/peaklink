function Threshold=getThreshold(XIC)
% We determin the threshold for the XIC automatically.
% First we determine noise variance.
XIC=XIC(XIC>0);
if isempty(XIC)==0
dnXIC=wden(XIC,'sqtwolog', 'h','mln',3,'db3');
mvavg=conv(XIC, ones(100,1)/100);
% cutout values greater then 500
mvavg=mvavg(mvavg<500);
[histmvavg, intx]=hist(mvavg,500);
[maxhist maxhistid]=max(histmvavg);
meanest=intx(maxhistid);
noisevar=var(XIC-dnXIC);
noisestd=sqrt(noisevar);

%maxint=max(XIC);

%if maxint/9 > 6*noisestd
%    Threshold=maxint/9;
%else
    Threshold=3*noisestd+meanest;
else
    Threshold=0;
end   
%end    