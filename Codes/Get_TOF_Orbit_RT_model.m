function [R_Training_ID01,T_TOF01_Orbit01,T_Training_ID01,R_TOF01_Orbit01,PP_TOF2Orbit01,PP_Normal_TOF2Orbit01,Mu_TOF2Orbit01,Sigma_TOF2Orbit01]=Get_TOF_Orbit_RT_model(Pep_com,TOF_ms2information_com,IntervalList_com01,retentiont01l1,score,Std_AMT_timeshift)

R_Training_ID01=[];
T_TOF01_Orbit01=[];
num=0;num1=0;
for i=1:length(Pep_com)
    Normalized_ms2time01_TOF=TOF_ms2information_com(i)/max(retentiont01l1);
    Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
    ID_T_inMS2interval=find(score{i}.Normal_T_TOF<=Normalized_ms2time01_TOF+3*Std_AMT_timeshift & score{i}.Normal_T_TOF>=Normalized_ms2time01_TOF-3*Std_AMT_timeshift);
    if ~isempty(ID_T_inMS2interval)
        num=num+1;
        V_AR=score{i}.R_TOFOrbit(ID_T_inMS2interval);
        AR_good_id=find(V_AR>=0.85);
        if ~isempty(AR_good_id)
%             if length(AR_good_id)==1
                T_T01=score{i}.T_TOF(ID_T_inMS2interval(AR_good_id));
                T_O01=score{i}.T_Orbit(ID_T_inMS2interval(AR_good_id));
                TD_T01O01=score{i}.T_TOFOrbit(ID_T_inMS2interval(AR_good_id));
                NT_T01=score{i}.Normal_T_TOF(ID_T_inMS2interval(AR_good_id));
                NT_O01=score{i}.Normal_T_Orbit(ID_T_inMS2interval(AR_good_id));
                NTD_T01O01=score{i}.Normal_T_TOFOrbit(ID_T_inMS2interval(AR_good_id));
                R_Training_ID01=[R_Training_ID01;i];
                T_TOF01_Orbit01=[T_TOF01_Orbit01; T_T01 T_O01 TD_T01O01 NT_T01 NT_O01 NTD_T01O01];
%             else num1=num1+1;
%             end
        end
    end
end

PP_TOF2Orbit01=polyfit(T_TOF01_Orbit01(:,1),T_TOF01_Orbit01(:,2),4);
PP_Normal_TOF2Orbit01=polyfit(T_TOF01_Orbit01(:,4),T_TOF01_Orbit01(:,5),4);
[Mu_TOF2Orbit01,Sigma_TOF2Orbit01]=normfit(polyval(PP_Normal_TOF2Orbit01,T_TOF01_Orbit01(:,4))-T_TOF01_Orbit01(:,5));
%%%%%%%%%%%%%%%%%% Get T_Training_data based on TOF01 and Orbit01
%%%%%%%%%%%%%%%%%%
T_Training_ID01=[];
R_TOF01_Orbit01=[];
num=0;num1=0;
for i=1:length(Pep_com)
    Normalized_ms2time01_TOF=TOF_ms2information_com(i)/max(retentiont01l1);
    Interval_matrix01=IntervalList_com01{i}.intervallist_after7_combine;
    ID_T_inMS2interval=find(score{i}.Normal_T_TOF<=Normalized_ms2time01_TOF+3*Std_AMT_timeshift & score{i}.Normal_T_TOF>=Normalized_ms2time01_TOF-3*Std_AMT_timeshift);
    if ~isempty(ID_T_inMS2interval)
        num=num+1;
        V_AT=polyval(PP_Normal_TOF2Orbit01,score{i}.Normal_T_TOF(ID_T_inMS2interval))-score{i}.Normal_T_Orbit(ID_T_inMS2interval);
        AT_good_id=find(abs(V_AT)<=Std_AMT_timeshift);
        if ~isempty(AT_good_id)
%             if length(AT_good_id)==1
                R_TOF01_Orbit01=[R_TOF01_Orbit01;score{i}.R_TOFOrbit(ID_T_inMS2interval(AT_good_id))];
                T_Training_ID01=[T_Training_ID01;i];
%             else num1=num1+1;
%             end            
        end
    end
end









