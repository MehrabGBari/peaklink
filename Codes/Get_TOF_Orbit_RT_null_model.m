function [T_null_TOF01_Orbit01,R_null_TOF01_Orbit01,PP_TOF2Orbit01,PP_Normal_TOF2Orbit01,Mu_TOF2Orbit01,Sigma_TOF2Orbit01,PHAT,T_null_95per_Bond,R_null_95per_Bond,Null_model_id]=Get_TOF_Orbit_RT_null_model(Pep_com,TOF_ms2information_com,IntervalList_com01,retentiont01l1,score,Std_AMT_timeshift)

R_null_TOF01_Orbit01=[];
T_null_TOF01_Orbit01=[];
LE_null_TOF01_Orbit01=[];
R_TOF01_Orbit01=[];
T_TOF01_Orbit01=[];
LE_TOF01_Orbit01=[];
Null_model_id=[];
num=0;num1=0;
for i=1:length(Pep_com)
    Normalized_ms2time01_TOF=TOF_ms2information_com(i)/max(retentiont01l1);
    Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
    ID_T_inMS2interval_1std=find(score{i}.Normal_T_TOF<=Normalized_ms2time01_TOF+1*Std_AMT_timeshift & score{i}.Normal_T_TOF>=Normalized_ms2time01_TOF-1*Std_AMT_timeshift);
%     ID_T_inMS2interval_3std=find(score{i}.Normal_T_TOF<=Normalized_ms2time01_TOF+3*Std_AMT_timeshift & score{i}.Normal_T_TOF>=Normalized_ms2time01_TOF-3*Std_AMT_timeshift);
    ID_T_inMS2interval_3std=1:length(score{i}.Normal_T_TOF);
    if ~isempty(ID_T_inMS2interval_1std) && length(ID_T_inMS2interval_1std)==1
        num=num+1;
        if length(ID_T_inMS2interval_3std)>=2
            num1=num1+1;
        end
        V_AR=score{i}.R_TOFOrbit(ID_T_inMS2interval_1std);
        if V_AR>=0.80 && length(ID_T_inMS2interval_3std)>=2
            
            Null_model_id=[Null_model_id;i];
            ID_T_null_inMS2interval_3std=ID_T_inMS2interval_3std;
            ID_T_null_inMS2interval_3std(ID_T_inMS2interval_3std==ID_T_inMS2interval_1std)=[];
            
            T_null_T01=score{i}.T_TOF(ID_T_null_inMS2interval_3std);
            T_null_O01=score{i}.T_Orbit(ID_T_null_inMS2interval_3std);
            TD_null_T01O01=score{i}.T_TOFOrbit(ID_T_null_inMS2interval_3std);
            NT_null_T01=score{i}.Normal_T_TOF(ID_T_null_inMS2interval_3std);
            NT_null_O01=score{i}.Normal_T_Orbit(ID_T_null_inMS2interval_3std);
            NTD_null_T01O01=score{i}.Normal_T_TOFOrbit(ID_T_null_inMS2interval_3std);
            LE_null_T01=score{i}.LE_TOF(ID_T_null_inMS2interval_3std);
            LE_null_T01O01=score{i}.LE_TOFOrbit(ID_T_null_inMS2interval_3std);            

            T_null_TOF01_Orbit01=[T_null_TOF01_Orbit01; T_null_T01 T_null_O01 TD_null_T01O01 NT_null_T01 NT_null_O01 NTD_null_T01O01];
            R_null_TOF01_Orbit01=[R_null_TOF01_Orbit01;score{i}.R_TOFOrbit(ID_T_null_inMS2interval_3std)];
            LE_null_TOF01_Orbit01=[LE_null_TOF01_Orbit01;score{i}.LE_TOFOrbit(ID_T_null_inMS2interval_3std)];
            
            T_T01=score{i}.T_TOF(ID_T_inMS2interval_1std);
            T_O01=score{i}.T_Orbit(ID_T_inMS2interval_1std);
            TD_T01O01=score{i}.T_TOFOrbit(ID_T_inMS2interval_1std);
            NT_T01=score{i}.Normal_T_TOF(ID_T_inMS2interval_1std);
            NT_O01=score{i}.Normal_T_Orbit(ID_T_inMS2interval_1std);
            NTD_T01O01=score{i}.Normal_T_TOFOrbit(ID_T_inMS2interval_1std);
            LE_T01=score{i}.LE_TOF(ID_T_inMS2interval_3std);
            LE_T01O01=score{i}.LE_TOFOrbit(ID_T_inMS2interval_3std);

            T_TOF01_Orbit01=[T_TOF01_Orbit01; T_T01 T_O01 TD_T01O01 NT_T01 NT_O01 NTD_T01O01];
            R_TOF01_Orbit01=[R_TOF01_Orbit01;score{i}.R_TOFOrbit(ID_T_inMS2interval_1std)];
            LE_TOF01_Orbit01=[LE_TOF01_Orbit01;score{i}.LE_TOFOrbit(ID_T_inMS2interval_1std)];
            
        end    
 
    end

end


length(Null_model_id)
PP_null_TOF2Orbit01=polyfit(T_null_TOF01_Orbit01(:,1),T_null_TOF01_Orbit01(:,2),4);
PP_null_Normal_TOF2Orbit01=polyfit(T_null_TOF01_Orbit01(:,4),T_null_TOF01_Orbit01(:,5),4);
[Mu_null_TOF2Orbit01,Sigma_null_TOF2Orbit01]=normfit(polyval(PP_null_Normal_TOF2Orbit01,T_null_TOF01_Orbit01(:,4))-T_null_TOF01_Orbit01(:,5));

Diff=1;
Up_B=1; Low_B=-1;
while abs(Diff)>=0.0001
    x=(Up_B+Low_B)/2;
    y=normcdf(x,Mu_TOF2Orbit01,Sigma_TOF2Orbit01);
    Diff=y-0.525;
    if Diff<=0
        Low_B=x;
    else
        Up_B=x;
    end
end
T_null_95per_Bond=[2*Mu_TOF2Orbit01-x,x];
% normcdf(0.0113,Mu_TOF2Orbit01,Sigma_TOF2Orbit01)

PHAT_null=gamfit(1-R_null_TOF01_Orbit01);

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






