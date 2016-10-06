function est=O18rateLinearPeak(p,intensity,min_f,max_f,rangerate,pp)
global p1;
global p2;
global p3;
global sigma1;
global sigma2;
global sigma3;
global yintensity;
global ns;
ns=0;
yintensity=intensity;
p1=p;
p2=zeros(1,6);p2(3:end)=p1(1:end-2);
p3=zeros(1,6);p3(5:end)=p1(1:end-4);
%p=[p1,p2,p3,p4,p5,p6];
sigma1=diag(p1)-p1'*p1;
sigma2=zeros(6);sigma2(3:end,3:end)=sigma1(1:end-2,1:end-2);
sigma3=zeros(6);sigma3(5:end,5:end)=sigma1(1:end-4,1:end-4);

init_a=(yintensity(1)*p1(1)+yintensity(2)*p1(2))/(p1(1)^2+p1(2)^2);
init_rate=(yintensity(5)-p1(5)/p1(1)*yintensity(1)-p1(3)/p1(1)*(yintensity(3)-p1(3)/p1(1)*yintensity(1))+(yintensity(3)-p1(3)/p1(1)*yintensity(1)))/yintensity(1);
init_f=1;
% init_b=init_a*init_rate;
% init_b2=init_b*init_f^2;
% init_b1=2*init_b*init_f*(1-init_f);
% init_b0=init_b*(1-init_f)^2;
% init_a1=init_a+init_b0;
% init_a2=init_b1;
% init_a3=init_b2;

lowerrate=rangerate*init_rate; upperrate=1/rangerate*init_rate;
if init_rate>999 || init_rate<0.001 || ~isfinite(init_rate)
    init_rate=1; lowerrate=1/100; upperrate=100;
end
if init_a<0.001
    init_a=0.001;
end

%options=optimset('Algorithm','active-set','Display','off');
%est=fmincon(@myLSEpeakfunc, [init_a, init_rate, init_f],[],[],[],[],[rangerate*init_a lowerrate min_f],[1/rangerate*init_a upperrate 1],[],options);
C=[p1' p2' p3'];d=yintensity';
peakindex=1:6;
for i=0:2:4
    if yintensity(i+1)>yintensity(i+2)
        peakindex(i+2)=0;
    else
        peakindex(i+1)=0;
    end
end
peakindex(peakindex==0)=[];
C=C(peakindex,:);d=d(peakindex);
P=[1, 1-pp, (1-pp)^2; 0, pp, 2*pp*(1-pp); 0, 0, pp^2];
C=C*P;
x=lsqnonneg(C,d);
f=2*x(3)/(2*x(3)+x(2));
theta2=x(2)/(2*f*(1-f));
theta1=x(1)-(1-f)^2*theta2;
%est=[theta1, theta2/theta1,f];
est=[theta1, theta2/theta1,f*pp];

if ~isfinite(est(2))
    est(2)=-1;
end
if ~isfinite(est(3))
    est(3)=-1;
end

if est(2)<0 || est(3)<min_f || est(3)>max_f
    options=optimset('Algorithm','active-set','Display','off');
    try
        est=fmincon(@myLSEpeakfunc, [init_a, init_rate,init_f],[],[],[],[],[rangerate*init_a lowerrate min_f],[1/rangerate*init_a upperrate max_f],[],options);
    catch exception
        rangerate=0.7;lowerrate=lowerrate*0.8;upperrate=upperrate*1.25;
        est=fmincon(@myLSEpeakfunc, [init_a, init_rate,init_f],[],[],[],[],[rangerate*init_a lowerrate min_f],[1/rangerate*init_a upperrate max_f],[],options);
        disp('one pass');
    end
end
if ~isfinite(est(2))
    est(2)=-1;
end
if ~isfinite(est(3))
    est(3)=-1;
end

end
