function [c, d, newc, newd, posi]=resampleforhalf(a, b)
% a=ElutionPeak01';
% b=ElutionPeak0102';

%%%%%%%%%
Rat_len=(length(b)-1)/(length(a)-1);
if Rat_len>=1
        c=b;   d=a;
        orignal=[zeros(1,length(a)),b,zeros(1,length(a))];
        filp=a;
        filpv1=filp(end:-1:1);
        C=conv(orignal,filpv1);
        [V,ID]=max(C);
%         newc=a; newd=orignal(ID-length(a)+1:ID);  %%% cut the longer one
        newc=[zeros(1,ID-length(a)),a,zeros(1,length(orignal)-ID)]; newd=orignal; %%% add zeros to the shorter one
        posi=0;
else
    c=a;    d=b;
    orignal=[zeros(1,length(b)),a,zeros(1,length(b))];
    filp=b;
    filpv1=filp(end:-1:1);
    C=conv(orignal,filpv1);
    [V,ID]=max(C);
%     newc=orignal(ID-length(b)+1:ID); newd=b; %%% cut the longer one
    newc=orignal; newd=[zeros(1,ID-length(b)),b,zeros(1,length(orignal)-ID)]; %%% add zeros to the shorter one
    posi=1;
end

%%%%%%%%%


%%%%%%%%%
% Rat_len=(length(b)-1)/(length(a)-1);
% if Rat_len>=1
%     YY = interp1([1:length(a)],a,[1:1/Rat_len:length(a)]);
%     c=b;   d=a;
%     orignal=YY-mean(YY);
%     filp=b-mean(b);
%     filpv1=filp(end:-1:1);
%     C=conv(orignal,filpv1);
%     [V,ID]=max(C);
%     if ID>=length(YY)
%         YYv1=[YY,zeros(1,ID-length(YY))];
%         bv1=[zeros(1,ID-length(YY)),b];
%     else
%         YYv1=[zeros(1,ID-length(YY)),YY];
%         bv1=[b,zeros(1,ID-length(YY))];
%     end
%     newc=YYv1; newd=bv1;
%     posi=0;
% else
%     YY = interp1([1:length(b)],b,[1:Rat_len:length(b)]);
%     c=a;    d=b;
%     orignal=a-mean(a);
%     filp=YY-mean(YY);
%     filpv1=filp(end:-1:1);
%     C=conv(orignal,filpv1);
%     [V,ID]=max(C);
%     if ID>=length(a)
%         av1=[a,zeros(1,ID-length(a))];
%         YYv1=[zeros(1,ID-length(a)),YY];
%     else
%         av1=[zeros(1,ID-length(a)),a];
%         YYv1=[YY,zeros(1,ID-length(a))];
%     end
%     
%     newc=av1; newd=YYv1;   
%     posi=1;
% end

%%%%%%%%%%%%%%%%%%%%%%

% am=max(a);
% bm=max(b);
% 
% aposi=find(a>=0.5*am);
% bposi=find(b>=0.5*bm);
% 
% ainterval=aposi(end)-aposi(1)+1;
% binterval=bposi(end)-bposi(1)+1;
% 
% vect=[ainterval,binterval];
% [N,posi]=max(vect);
% 
% if posi==1
% c=a;    d=b;   cposi=aposi;  dposi=bposi; e=2;
% else c=b;   d=a;   cposi=bposi;  dposi=aposi; e=1;
% end
% 
% d1=pchip([1:length(d)],d,[1:vect(e)/N:length(d)]);
% 
% max_pc=cposi(1);
% 
% d1m=max(d1);
% d1posi=find(d1>=0.5*d1m);
% d1interval=d1posi(end)-d1posi(1)+1;
% max_pd1=d1posi(1);
% 
% front=max(max_pc,max_pd1);
% back=max(length(c)-max_pc+1, length(d1)-max_pd1+1);
% newc=zeros(1,back+front);
% newd=zeros(1,back+front);
% 
% newc(front-max_pc+1:length(c)+front-max_pc)=c;
% newd(front-max_pd1+1:length(d1)+front-max_pd1)=d1;





