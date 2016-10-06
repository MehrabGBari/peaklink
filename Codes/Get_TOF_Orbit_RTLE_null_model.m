function [T_null_TOF_Orbit,R_null_TOF_Orbit,LE_null_TOF_Orbit...
                PP_null_TOF2Orbit,PP_null_Normal_TOF2Orbit,...
                Mu_null_TOF2Orbit,Sigma_null_TOF2Orbit,...
                PHAT_null,T_null_95per_Bond,R_null_95per_Bond]=Get_TOF_Orbit_RTLE_null_model(Training_Orbitms2time,...
                        Training_Intervallist,retentiontl1,score,Std_AMT_timeshift,Training_TOFInterval_id)
                    
R_null_TOF_Orbit=[];
T_null_TOF_Orbit=[];
LE_null_TOF_Orbit=[];
R_TOF_Orbit=[];
T_TOF_Orbit=[];
LE_TOF_Orbit=[];
Null_model_id=[];
num=0;num1=0;
for i=1:length(Training_Orbitms2time)
    Training_Intervallist{i}.Labelling_efficiency_after7_combine;
    k=Training_TOFInterval_id(i);
    if score{i}.T_TOF(1)~=0
        num=num+1;
        T_T=score{i}.T_TOF(k);
        T_O=score{i}.T_Orbit(k);
        TD_TO=score{i}.T_TOFOrbit(k);
        NT_T=score{i}.Normal_T_TOF(k);
        NT_O=score{i}.Normal_T_Orbit(k);
        NTD_TO=score{i}.Normal_T_TOFOrbit(k);
        LE_T=score{i}.LE_TOF(k);
        LE_TO=score{i}.LE_TOFOrbit(k);

        T_TOF_Orbit=[T_TOF_Orbit; T_T T_O TD_TO NT_T NT_O NTD_TO];
        R_TOF_Orbit=[R_TOF_Orbit;score{i}.R_TOFOrbit(k)];
        LE_TOF_Orbit=[LE_TOF_Orbit;score{i}.LE_TOFOrbit(k)];

        ID=1:length(score{i}.Normal_T_TOF);
        ID(k)=[];
        if ~isempty(ID)
                num1=num1+1;
                T_null_T=score{i}.T_TOF(ID);
                T_null_O=score{i}.T_Orbit(ID);
                TD_null_TO=score{i}.T_TOFOrbit(ID);
                NT_null_T=score{i}.Normal_T_TOF(ID);
                NT_null_O=score{i}.Normal_T_Orbit(ID);
                NTD_null_TO=score{i}.Normal_T_TOFOrbit(ID);
                LE_null_T=score{i}.LE_TOF(ID);
                LE_null_TO=score{i}.LE_TOFOrbit(ID);            

                T_null_TOF_Orbit=[T_null_TOF_Orbit; T_null_T T_null_O TD_null_TO NT_null_T NT_null_O NTD_null_TO];
                R_null_TOF_Orbit=[R_null_TOF_Orbit;score{i}.R_TOFOrbit(ID)];
                LE_null_TOF_Orbit=[LE_null_TOF_Orbit;score{i}.LE_TOFOrbit(ID)'];
        end
    end
    
%     Normalized_ms2time01_TOF=Training_Orbitms2time01(i)/max(retentiont01l1);
%     Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
%     ID_T_inMS2interval_1std=find(score{i}.Normal_T_TOF<=Normalized_ms2time01_TOF+1*Std_AMT_timeshift & score{i}.Normal_T_TOF>=Normalized_ms2time01_TOF-1*Std_AMT_timeshift);
% %     ID_T_inMS2interval_3std=find(score{i}.Normal_T_TOF<=Normalized_ms2time01_TOF+3*Std_AMT_timeshift & score{i}.Normal_T_TOF>=Normalized_ms2time01_TOF-3*Std_AMT_timeshift);
%     ID_T_inMS2interval_3std=1:length(score{i}.Normal_T_TOF);
%     if ~isempty(ID_T_inMS2interval_1std) && length(ID_T_inMS2interval_1std)==1
%         num=num+1;
%         if length(ID_T_inMS2interval_3std)>=2
%             num1=num1+1;
%         end
%         V_AR=score{i}.R_TOFOrbit(ID_T_inMS2interval_1std);
%         if V_AR>=0.80 && length(ID_T_inMS2interval_3std)>=2
%             
%             Null_model_id=[Null_model_id;i];
%             ID_T_null_inMS2interval_3std=ID_T_inMS2interval_3std;
%             ID_T_null_inMS2interval_3std(ID_T_inMS2interval_3std==ID_T_inMS2interval_1std)=[];
%             
%             T_null_T01=score{i}.T_TOF(ID_T_null_inMS2interval_3std);
%             T_null_O01=score{i}.T_Orbit(ID_T_null_inMS2interval_3std);
%             TD_null_T01O01=score{i}.T_TOFOrbit(ID_T_null_inMS2interval_3std);
%             NT_null_T01=score{i}.Normal_T_TOF(ID_T_null_inMS2interval_3std);
%             NT_null_O01=score{i}.Normal_T_Orbit(ID_T_null_inMS2interval_3std);
%             NTD_null_T01O01=score{i}.Normal_T_TOFOrbit(ID_T_null_inMS2interval_3std);
%             LE_null_T01=score{i}.LE_TOF(ID_T_null_inMS2interval_3std);
%             LE_null_T01O01=score{i}.LE_TOFOrbit(ID_T_null_inMS2interval_3std);            
% 
%             T_null_TOF01_Orbit01=[T_null_TOF01_Orbit01; T_null_T01 T_null_O01 TD_null_T01O01 NT_null_T01 NT_null_O01 NTD_null_T01O01];
%             R_null_TOF01_Orbit01=[R_null_TOF01_Orbit01;score{i}.R_TOFOrbit(ID_T_null_inMS2interval_3std)];
%             LE_null_TOF01_Orbit01=[LE_null_TOF01_Orbit01;score{i}.LE_TOFOrbit(ID_T_null_inMS2interval_3std)];
%             
%             T_T01=score{i}.T_TOF(ID_T_inMS2interval_1std);
%             T_O01=score{i}.T_Orbit(ID_T_inMS2interval_1std);
%             TD_T01O01=score{i}.T_TOFOrbit(ID_T_inMS2interval_1std);
%             NT_T01=score{i}.Normal_T_TOF(ID_T_inMS2interval_1std);
%             NT_O01=score{i}.Normal_T_Orbit(ID_T_inMS2interval_1std);
%             NTD_T01O01=score{i}.Normal_T_TOFOrbit(ID_T_inMS2interval_1std);
%             LE_T01=score{i}.LE_TOF(ID_T_inMS2interval_3std);
%             LE_T01O01=score{i}.LE_TOFOrbit(ID_T_inMS2interval_3std);
% 
%             T_TOF01_Orbit01=[T_TOF01_Orbit01; T_T01 T_O01 TD_T01O01 NT_T01 NT_O01 NTD_T01O01];
%             R_TOF01_Orbit01=[R_TOF01_Orbit01;score{i}.R_TOFOrbit(ID_T_inMS2interval_1std)];
%             LE_TOF01_Orbit01=[LE_TOF01_Orbit01;score{i}.LE_TOFOrbit(ID_T_inMS2interval_1std)];
%             
%         end    
%  
%     end

end

figure
hist(R_null_TOF_Orbit)
figure
hist(R_TOF_Orbit)
figure
hist(T_null_TOF_Orbit(:,4)-T_null_TOF_Orbit(:,5))
figure
hist(T_TOF_Orbit(:,4)-T_TOF_Orbit(:,5))
figure
hist(abs(LE_null_TOF_Orbit),25)
figure
hist(abs(LE_TOF_Orbit),25)

PP_null_TOF2Orbit=polyfit(T_null_TOF_Orbit(:,1),T_null_TOF_Orbit(:,2),4);
PP_null_Normal_TOF2Orbit=polyfit(T_null_TOF_Orbit(:,4),T_null_TOF_Orbit(:,5),4);
[Mu_null_TOF2Orbit,Sigma_null_TOF2Orbit]=normfit(polyval(PP_null_Normal_TOF2Orbit,T_null_TOF_Orbit(:,4))-T_null_TOF_Orbit(:,5));

Diff=1;
Up_B=1; Low_B=-1;
while abs(Diff)>=0.0001
    x=(Up_B+Low_B)/2;
    y=normcdf(x,Mu_TOF2Orbit,Sigma_TOF2Orbit);
    Diff=y-0.525;
    if Diff<=0
        Low_B=x;
    else
        Up_B=x;
    end
end
T_null_95per_Bond=[2*Mu_TOF2Orbit-x,x];
% normcdf(0.0113,Mu_TOF2Orbit01,Sigma_TOF2Orbit01)

PHAT_null=gamfit(1-R_null_TOF_Orbit);

Diff=1;
Up_B=1; Low_B=0;
while abs(Diff)>=0.0001
    x=(Up_B+Low_B)/2;
    y=gamcdf(x,PHAT_null(1),PHAT_null(2));
    Diff=y-0.05;
    if Diff>=0
        Up_B=x;
    else
        Low_B=x;
    end
end

R_null_95per_Bond=[1-x,1];






