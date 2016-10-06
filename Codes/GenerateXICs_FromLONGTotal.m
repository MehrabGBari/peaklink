function TOF_XICs=GenerateXICs_FromLONGTotal(TOF_XICs,XICs_C_N,XICs_R_N,TOFAMT_LIST,FractionName,dataname)


Fract=find(TOFAMT_LIST(:,3)==FractionName);
% Fract_cs1=find(TOFAMT_LIST(Fract,2)==1);
Fract_cs2=find(TOFAMT_LIST(Fract,2)==2);
Fract_cs3=find(TOFAMT_LIST(Fract,2)==3);
Fract_cs4=find(TOFAMT_LIST(Fract,2)==4);
Fract_cs5=find(TOFAMT_LIST(Fract,2)==5);

CS=2;
TOF_XICs=GenerateXICs_atdifferentcharge(TOF_XICs,XICs_R_N,Fract,Fract_cs2,TOFAMT_LIST,FractionName,dataname,CS);
CS=3;
TOF_XICs=GenerateXICs_atdifferentcharge(TOF_XICs,XICs_R_N,Fract,Fract_cs3,TOFAMT_LIST,FractionName,dataname,CS);
CS=4;
TOF_XICs=GenerateXICs_atdifferentcharge(TOF_XICs,XICs_R_N,Fract,Fract_cs4,TOFAMT_LIST,FractionName,dataname,CS);
CS=5;
TOF_XICs=GenerateXICs_atdifferentcharge(TOF_XICs,XICs_R_N,Fract,Fract_cs5,TOFAMT_LIST,FractionName,dataname,CS);

% if ~isempty(Fract_cs2)
%     C_ID=TOFAMT_LIST(Fract(Fract_cs2),4);
%     filename=['D:\Program\QTOF_replicate_identification\LongXIC5574f4\Allinteract', num2str(FractionName) ,'_3.ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_0', num2str(dataname),'_f4.centroid1.tolerance20.usemax0.chargeXIC2.mat'];
%     load(filename);
%     for i=1:6
%         TOF_XICs{i}(:, C_ID)=chargeXIC2{i}(1:XICs_R_N, C_ID);
%     end
% end
% if ~isempty(Fract_cs3)
%     C_ID=TOFAMT_LIST(Fract(Fract_cs3),4);
%     filename=['D:\Program\QTOF_replicate_identification\LongXIC5574f4\Allinteract', num2str(FractionName) ,'_3.ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_0', num2str(dataname),'_f4.centroid1.tolerance20.usemax0.chargeXIC3.mat'];
%     load(filename);
%     for i=1:6
%         TOF_XICs{i}(:, C_ID)=chargeXIC3{i}(1:XICs_R_N, C_ID);
%     end
% end
% if ~isempty(Fract_cs4)
%     C_ID=TOFAMT_LIST(Fract(Fract_cs4),4);
%     filename=['D:\Program\QTOF_replicate_identification\LongXIC5574f4\Allinteract', num2str(FractionName) ,'_3.ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_0', num2str(dataname),'_f4.centroid1.tolerance20.usemax0.chargeXIC4.mat'];
%     load(filename);
%     for i=1:6
%         TOF_XICs{i}(:, C_ID)=chargeXIC4{i}(1:XICs_R_N, C_ID);
%     end
% end
% if ~isempty(Fract_cs5)
%     C_ID=TOFAMT_LIST(Fract(Fract_cs5),4);
%     filename=['D:\Program\QTOF_replicate_identification\LongXIC5574f4\Allinteract', num2str(FractionName) ,'_3.ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_0', num2str(dataname),'_f4.centroid1.tolerance20.usemax0.chargeXIC5.mat'];
%     load(filename);
%     for i=1:6
%         TOF_XICs{i}(:, C_ID)=chargeXIC5{i}(1:XICs_R_N, C_ID);
%     end
% end










