function [PP,T_mu_corr,T_sigma_corr,T_mu_non_corr,T_sigma_non_corr,...
                LE_pdf_pp_corr,LE_pdf_pp_non_corr,...
                KL_pdf_pp_corr,KL_pdf_pp_non_corr,...
                LEKLT_pdf_pp_corr,LEKLT_pdf_pp_non_corr,...
                LET_pdf_pp_corr,LET_pdf_pp_non_corr,...
                LEKL_pdf_pp_corr,LEKL_pdf_pp_non_corr,...
                KLT_pdf_pp_corr,KLT_pdf_pp_non_corr,...
                LE_cdf_pp_corr,LE_cdf_pp_non_corr,...
                KL_cdf_pp_corr,KL_cdf_pp_non_corr,...
                LEKLT_cdf_pp_corr,LEKLT_cdf_pp_non_corr,...
                LET_cdf_pp_corr,LET_cdf_pp_non_corr,...
                LEKL_cdf_pp_corr,LEKL_cdf_pp_non_corr,...
                KLT_cdf_pp_corr,KLT_cdf_pp_non_corr]=Generate_RTLEmodel_Parameters(R_corr_TOF_Orbit,R_non_corr_TOF_Orbit,...
                    T_corr_TOF_Orbit,T_non_corr_TOF_Orbit,...
                    LE_corr_TOF_Orbit,LE_non_corr_TOF_Orbit,...
                    KL_corr_TOF_Orbit,KL_non_corr_TOF_Orbit)

% ID_R_good=find(R_corr_TOF_Orbit>=0.90);
% LE_corr_TOF_Orbitv1=LE_corr_TOF_Orbit(ID_R_good);
% ID_LE_good=find(abs(LE_corr_TOF_Orbit)<=0.05);
% R_corr_TOF_Orbitv1=R_corr_TOF_Orbit(ID_LE_good);

PP=polyfit(T_corr_TOF_Orbit(:,1),T_corr_TOF_Orbit(:,2),4);
T_corr_diff=polyval(PP,T_corr_TOF_Orbit(:,1))-T_corr_TOF_Orbit(:,2);
T_non_corr_diff=polyval(PP,T_non_corr_TOF_Orbit(:,1))-T_non_corr_TOF_Orbit(:,2);
xxx=0:0.01:1;
figure
plot(T_corr_TOF_Orbit(:,1),T_corr_TOF_Orbit(:,2),'r*')
hold on
plot(xxx,polyval(PP,xxx))
grid on
xlabel('Normalized Time in Data 01');
ylabel('Normalized Time in Data 02');
legend('Time point','Gwarping function');

[T_mu_corr,T_sigma_corr]=normfit(T_corr_diff);
[T_mu_non_corr,T_sigma_non_corr]=normfit(T_non_corr_diff);
% T_corr_PARMHAT=gamfit(abs(T_corr_diff));
% T_non_corr_PARMHAT=gamfit(abs(T_non_corr_diff));

x=-1:0.01:1;
figure
subplot(2,2,1)
hist(T_corr_diff);grid on;title('(a)');xlabel('Normalized Time shift');
subplot(2,2,3)
hist(T_non_corr_diff);grid on;title('(b)');xlabel('Normalized Time shift');
subplot(1,2,2)
plot(x,normpdf(x,T_mu_corr,T_sigma_corr),'r');
hold on
plot(x,normpdf(x,T_mu_non_corr,T_sigma_non_corr));
grid on;title('(c)');xlabel('Normalized Time shift');
legend('Corresponding Time shift','Non-Corr Tme shift')

LE_corr_TOF_Orbitv1=LE_corr_TOF_Orbit(~isinf(LE_corr_TOF_Orbit));
LE_non_corr_TOF_Orbitv1=LE_non_corr_TOF_Orbit(~isinf(LE_non_corr_TOF_Orbit));
% LE_corr_PARMHAT=gamfit(abs(LE_corr_TOF_Orbitv1));
% LE_non_corr_PARMHAT=gamfit(abs(LE_non_corr_TOF_Orbitv1));

X_center=-3:0.2:4;
[Y_corr,I_corr]=hist((LE_corr_TOF_Orbitv1),X_center);
[Y_non_corr,I_non_corr]=hist((LE_non_corr_TOF_Orbitv1),X_center);

% Y_w_corr=wden(10*Y_corr./sum(Y_corr),'sqtwolog','h','mln',8,'sym8');
% Y_w_non_corr=wden(10*Y_non_corr./sum(Y_non_corr),'sqtwolog','h','mln',8,'sym8');
Y_w_corr=5*Y_corr./sum(Y_corr);
Y_w_non_corr=5*Y_non_corr./sum(Y_non_corr);

LE_pdf_pp_corr=pchip(I_corr,Y_w_corr);
LE_pdf_pp_non_corr=pchip(I_non_corr,Y_w_non_corr);

for i=1:length(Y_corr)
    LE_Number_corr_cdffit(i)=sum(Y_corr(1:i))/sum(Y_corr);
end
for i=1:length(Y_non_corr)
    LE_Number_non_corr_cdffit(i)=sum(Y_non_corr(1:i))/sum(Y_non_corr);
end
LE_cdf_pp_corr=pchip(I_corr,LE_Number_corr_cdffit);
LE_cdf_pp_non_corr=pchip(I_non_corr,LE_Number_non_corr_cdffit);
figure
subplot(2,2,1)
hist(LE_corr_TOF_Orbitv1);grid on;title('(a)');xlabel('Labeling difference');
subplot(2,2,3)
hist(LE_non_corr_TOF_Orbitv1);grid on;title('(b)');xlabel('Labeling difference');
subplot(1,2,2)
plot(I_corr,ppval(LE_pdf_pp_corr,I_corr),'r*-')
hold on
plot(I_non_corr,ppval(LE_pdf_pp_non_corr,I_non_corr),'k*-')
grid on;title('(c)');legend('Corresponding LE difference','Non-corr LE difference');
xlabel('Labeling difference');

% subplot(2,2,2)
% plot(X_center,10*Y_corr./sum(Y_corr),'ko-')
% hold on
% plot(I_corr,ppval(LE_pdf_pp_corr,I_corr),'r*-')
% grid on
% subplot(2,2,4)
% plot(X_center,10*Y_non_corr./sum(Y_non_corr),'go-')
% hold on
% plot(I_non_corr,ppval(LE_pdf_pp_non_corr,I_non_corr),'b*-')
% grid on

% R_corr_PARMHAT=betafit(R_corr_TOF_Orbit);
% R_non_corr_PARMHAT=betafit(R_non_corr_TOF_Orbit);
% 
% x=0:0.01:1;
% figure
% subplot(2,2,1)
% hist(R_corr_TOF_Orbit,50)
% subplot(2,2,3)
% hist(R_non_corr_TOF_Orbit,50)
% subplot(1,2,2)
% plot(x,betapdf(x,R_corr_PARMHAT(1),R_corr_PARMHAT(1)),'r');
% hold on
% plot(x,betapdf(x,R_non_corr_PARMHAT(1),R_non_corr_PARMHAT(1)));
% grid on


X_center=-25:1:0;
[Y_corr,I_corr]=hist(KL_corr_TOF_Orbit,X_center);
[Y_non_corr,I_non_corr]=hist(KL_non_corr_TOF_Orbit,X_center);

% Y_w_corr=wden(10*Y_corr./sum(Y_corr),'sqtwolog','h','sln',1,'sym8');
% Y_w_non_corr=wden(10*Y_non_corr./sum(Y_non_corr),'sqtwolog','h','sln',1,'sym8');
Y_w_corr=10*Y_corr./sum(Y_corr);
Y_w_non_corr=10*Y_non_corr./sum(Y_non_corr);

KL_pdf_pp_corr=pchip(I_corr,Y_w_corr);
KL_pdf_pp_non_corr=pchip(I_non_corr,Y_w_non_corr);
figure
subplot(2,1,1)
plot(X_center,10*Y_corr./sum(Y_corr),'ko-')
hold on
plot(I_corr,ppval(KL_pdf_pp_corr,I_corr),'r*-')
grid on
subplot(2,1,2)
plot(X_center,10*Y_non_corr./sum(Y_non_corr),'ko-')
hold on
plot(I_non_corr,ppval(KL_pdf_pp_non_corr,I_non_corr),'r*-')
grid on

figure
subplot(2,2,1)
hist(KL_corr_TOF_Orbit);grid on;title('(a)');xlabel('KL difference');
subplot(2,2,3)
hist(KL_non_corr_TOF_Orbit);grid on;title('(b)');xlabel('KL difference');
subplot(1,2,2)
plot(I_corr,ppval(KL_pdf_pp_corr,I_corr),'r*-')
hold on
plot(I_non_corr,ppval(KL_pdf_pp_non_corr,I_non_corr),'ko-')
grid on;
legend('Corresponding KL distance','Non-Corr KL distance');title('KL distance');



for i=1:length(Y_corr)
    KL_Number_corr_cdffit(i)=sum(Y_corr(1:i))/sum(Y_corr);
end
for i=1:length(Y_non_corr)
    KL_Number_non_corr_cdffit(i)=sum(Y_non_corr(1:i))/sum(Y_non_corr);
end
KL_cdf_pp_corr=pchip(I_corr,KL_Number_corr_cdffit);
KL_cdf_pp_non_corr=pchip(I_non_corr,KL_Number_non_corr_cdffit);
% figure
% plot(I_corr,ppval(pp_corr,I_corr),'r')
% hold on
% plot(I_non_corr,ppval(pp_non_corr,I_non_corr))


% KL_corr_PARMHAT=gamfit(-(KL_corr_TOF_Orbit));
% KL_non_corr_PARMHAT=gamfit(-(KL_non_corr_TOF_Orbit));
% [KL_corr_mu,KL_corr_sigma]=normfit((KL_corr_TOF_Orbit));
% [KL_non_corr_mu,KL_non_corr_sigma]=normfit((KL_non_corr_TOF_Orbit));
% 
% X_center=-20:1:0;
% [Y_corr,I_corr]=hist((KL_corr_TOF_Orbit),X_center);
% [Y_non_corr,I_non_corr]=hist((KL_non_corr_TOF_Orbit),X_center);
% figure
% subplot(2,2,1)
% hist((KL_corr_TOF_Orbit),X_center)
% grid on
% subplot(2,2,3)
% hist((KL_non_corr_TOF_Orbit),X_center)
% grid on
% subplot(2,2,2)
% plot(I_corr,Y_corr/sum(Y_corr),'k')
% hold on
% plot(X_center,normpdf(X_center,KL_corr_mu,KL_corr_sigma)./sum(normpdf(X_center,KL_corr_mu,KL_corr_sigma)),'r')
% grid on
% subplot(2,2,4)
% plot(I_non_corr,Y_non_corr/sum(Y_non_corr),'k')
% hold on
% plot(X_center,normpdf(X_center,KL_non_corr_mu,KL_non_corr_sigma)./sum(normpdf(X_center,KL_non_corr_mu,KL_non_corr_sigma)))
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Total KL T LE model
Total_LE_null_data=LE_non_corr_TOF_Orbit(~isinf(LE_non_corr_TOF_Orbit));
Total_T_null_data=T_non_corr_diff(~isinf(LE_non_corr_TOF_Orbit));
Total_KL_null_data=KL_non_corr_TOF_Orbit(~isinf(LE_non_corr_TOF_Orbit));
Total_LE_data=LE_corr_TOF_Orbit(~isinf(LE_corr_TOF_Orbit));
Total_T_data=T_corr_diff(~isinf(LE_corr_TOF_Orbit));
Total_KL_data=KL_corr_TOF_Orbit(~isinf(LE_corr_TOF_Orbit));
%%%%%%%%%%%%P means probility
P_total_LE_null_data=ppval(LE_pdf_pp_non_corr,Total_LE_null_data);
P_total_T_null_data=normpdf(Total_T_null_data,T_mu_non_corr,T_sigma_non_corr);
P_total_KL_null_data=ppval(KL_pdf_pp_non_corr,Total_KL_null_data);
P_total_LE_data=ppval(LE_pdf_pp_corr,Total_LE_data);
P_total_T_data=normpdf(Total_T_data,T_mu_corr,T_sigma_corr);
P_total_KL_data=ppval(KL_pdf_pp_corr,Total_KL_data);


%%%%%%%%%%%%%%%%% total LE T KL
X_cen=-12:1:12;
[Y_non_corr,I_non_corr]=hist(log(P_total_LE_null_data.*P_total_T_null_data.*P_total_KL_null_data),X_cen);
[Y_corr,I_corr]=hist(log(P_total_LE_data.*P_total_T_data.*P_total_KL_data),X_cen);
% Y_w_corr=wden(10*Y_corr./sum(Y_corr),'sqtwolog','h','mln',8,'sym8');
% Y_w_non_corr=wden(10*Y_non_corr./sum(Y_non_corr),'sqtwolog','h','mln',8,'sym8');
Y_w_corr=1*Y_corr./sum(Y_corr);
Y_w_non_corr=1*Y_non_corr./sum(Y_non_corr);
LEKLT_pdf_pp_corr=pchip(I_corr,Y_w_corr);
LEKLT_pdf_pp_non_corr=pchip(I_non_corr,Y_w_non_corr);
for i=1:length(Y_corr)
    LEKLT_Number_corr_cdffit(i)=sum(Y_corr(1:i))/sum(Y_corr);
end
for i=1:length(Y_non_corr)
    LEKLT_Number_non_corr_cdffit(i)=sum(Y_non_corr(1:i))/sum(Y_non_corr);
end
LEKLT_cdf_pp_corr=pchip(I_corr,LEKLT_Number_corr_cdffit);
LEKLT_cdf_pp_non_corr=pchip(I_non_corr,LEKLT_Number_non_corr_cdffit);

figure
plot(X_cen,ppval(LEKLT_cdf_pp_corr,X_cen),'r')
hold on
plot(X_cen,ppval(LEKLT_cdf_pp_non_corr,X_cen))
grid on;title('CDF');legend('Corresponding','Non-corresponding');

figure
subplot(2,2,1)
hist(log(P_total_LE_data.*P_total_T_data.*P_total_KL_data),X_cen) %%%%% LE T KL
grid on;title('(a)');
subplot(2,2,3)
hist(log(P_total_LE_null_data.*P_total_T_null_data.*P_total_KL_null_data),X_cen) %%%%% LE T KL
grid on;title('(b)');
subplot(1,2,2)
plot(I_corr,ppval(LEKLT_pdf_pp_corr,I_corr),'r*-');hold on;
plot(I_non_corr,ppval(LEKLT_pdf_pp_non_corr,I_non_corr),'k*-');grid on;
title('(c)');legend('Corresponding','Non-corresponding');

subplot(2,2,2)
plot(X_cen,1*Y_corr./sum(Y_corr),'ko')
hold on
plot(I_corr,ppval(LEKLT_pdf_pp_corr,I_corr),'r*-')
grid on
subplot(2,2,4)
plot(X_cen,1*Y_non_corr./sum(Y_non_corr),'ko')
hold on
plot(I_non_corr,ppval(LEKLT_pdf_pp_non_corr,I_non_corr),'r*-')
grid on

%%%%%%%%%%%%%%%%% LE T
X_cen=-12:1:12;
[Y_non_corr,I_non_corr]=hist(log(P_total_LE_null_data.*P_total_T_null_data),X_cen);
[Y_corr,I_corr]=hist(log(P_total_LE_data.*P_total_T_data),X_cen);
% Y_w_corr=wden(10*Y_corr./sum(Y_corr),'sqtwolog','h','mln',8,'sym8');
% Y_w_non_corr=wden(10*Y_non_corr./sum(Y_non_corr),'sqtwolog','h','mln',8,'sym8');
Y_w_corr=1*Y_corr./sum(Y_corr);
Y_w_non_corr=1*Y_non_corr./sum(Y_non_corr);
LET_pdf_pp_corr=pchip(I_corr,Y_w_corr);
LET_pdf_pp_non_corr=pchip(I_non_corr,Y_w_non_corr);
for i=1:length(Y_corr)
    LET_Number_corr_cdffit(i)=sum(Y_corr(1:i))/sum(Y_corr);
end
for i=1:length(Y_non_corr)
    LET_Number_non_corr_cdffit(i)=sum(Y_non_corr(1:i))/sum(Y_non_corr);
end
LET_cdf_pp_corr=pchip(I_corr,LET_Number_corr_cdffit);
LET_cdf_pp_non_corr=pchip(I_non_corr,LET_Number_non_corr_cdffit);

figure
subplot(2,2,1)
hist(log(P_total_LE_data.*P_total_T_data),X_cen) %%%%% LE T KL
grid on
subplot(2,2,3)
hist(log(P_total_LE_null_data.*P_total_T_null_data),X_cen) %%%%% LE T KL
grid on
subplot(2,2,2)
plot(X_cen,1*Y_corr./sum(Y_corr),'ko')
hold on
plot(I_corr,ppval(LET_pdf_pp_corr,I_corr),'r*-')
grid on
subplot(2,2,4)
plot(X_cen,1*Y_non_corr./sum(Y_non_corr),'ko')
hold on
plot(I_non_corr,ppval(LET_pdf_pp_non_corr,I_non_corr),'r*-')
grid on

%%%%%%%%%%%%%%%%% LE KL
X_cen=-12:1:12;
[Y_non_corr,I_non_corr]=hist(log(P_total_LE_null_data.*P_total_KL_null_data),X_cen);
[Y_corr,I_corr]=hist(log(P_total_LE_data.*P_total_KL_data),X_cen);
% Y_w_corr=wden(10*Y_corr./sum(Y_corr),'sqtwolog','h','mln',8,'sym8');
% Y_w_non_corr=wden(10*Y_non_corr./sum(Y_non_corr),'sqtwolog','h','mln',8,'sym8');
Y_w_corr=1*Y_corr./sum(Y_corr);
Y_w_non_corr=1*Y_non_corr./sum(Y_non_corr);
LEKL_pdf_pp_corr=pchip(I_corr,Y_w_corr);
LEKL_pdf_pp_non_corr=pchip(I_non_corr,Y_w_non_corr);
for i=1:length(Y_corr)
    LEKL_Number_corr_cdffit(i)=sum(Y_corr(1:i))/sum(Y_corr);
end
for i=1:length(Y_non_corr)
    LEKL_Number_non_corr_cdffit(i)=sum(Y_non_corr(1:i))/sum(Y_non_corr);
end
LEKL_cdf_pp_corr=pchip(I_corr,LEKL_Number_corr_cdffit);
LEKL_cdf_pp_non_corr=pchip(I_non_corr,LEKL_Number_non_corr_cdffit);

figure
subplot(2,2,1)
hist(log(P_total_LE_data.*P_total_KL_data),X_cen) %%%%% LE T KL
grid on
subplot(2,2,3)
hist(log(P_total_LE_null_data.*P_total_KL_null_data),X_cen) %%%%% LE T KL
grid on
subplot(2,2,2)
plot(X_cen,1*Y_corr./sum(Y_corr),'ko')
hold on
plot(I_corr,ppval(LEKL_pdf_pp_corr,I_corr),'r*-')
grid on
subplot(2,2,4)
plot(X_cen,1*Y_non_corr./sum(Y_non_corr),'ko')
hold on
plot(I_non_corr,ppval(LEKL_pdf_pp_non_corr,I_non_corr),'r*-')
grid on

%%%%%%%%%%%%%%%%% KL T
X_cen=-12:1:12;
[Y_non_corr,I_non_corr]=hist(log(P_total_T_null_data.*P_total_KL_null_data),X_cen);
[Y_corr,I_corr]=hist(log(P_total_T_data.*P_total_KL_data),X_cen);
% Y_w_corr=wden(10*Y_corr./sum(Y_corr),'sqtwolog','h','mln',8,'sym8');
% Y_w_non_corr=wden(10*Y_non_corr./sum(Y_non_corr),'sqtwolog','h','mln',8,'sym8');
Y_w_corr=1*Y_corr./sum(Y_corr);
Y_w_non_corr=1*Y_non_corr./sum(Y_non_corr);
KLT_pdf_pp_corr=pchip(I_corr,Y_w_corr);
KLT_pdf_pp_non_corr=pchip(I_non_corr,Y_w_non_corr);
for i=1:length(Y_corr)
    KLT_Number_corr_cdffit(i)=sum(Y_corr(1:i))/sum(Y_corr);
end
for i=1:length(Y_non_corr)
    KLT_Number_non_corr_cdffit(i)=sum(Y_non_corr(1:i))/sum(Y_non_corr);
end
KLT_cdf_pp_corr=pchip(I_corr,KLT_Number_corr_cdffit);
KLT_cdf_pp_non_corr=pchip(I_non_corr,KLT_Number_non_corr_cdffit);

figure
subplot(2,2,1)
hist(log(P_total_T_data.*P_total_KL_data),X_cen) %%%%% LE T KL
grid on
subplot(2,2,3)
hist(log(P_total_T_null_data.*P_total_KL_null_data),X_cen) %%%%% LE T KL
grid on
subplot(2,2,2)
plot(X_cen,1*Y_corr./sum(Y_corr),'ko')
hold on
plot(I_corr,ppval(KLT_pdf_pp_corr,I_corr),'r*-')
grid on
subplot(2,2,4)
plot(X_cen,1*Y_non_corr./sum(Y_non_corr),'ko')
hold on
plot(I_non_corr,ppval(KLT_pdf_pp_non_corr,I_non_corr),'r*-')
grid on

% %%%%%%%%%%%%%%%%%
% figure
% subplot(2,1,1)
% hist(log(P_total_T_null_data),X_cen) %%%%% LE T KL
% grid on
% subplot(2,1,2)
% hist(log(P_total_T_data),X_cen) %%%%% LE T KL
% grid on
% %%%%%%%%%%%%%%%%%
