function TOF_XICs=GenerateXICs_atdifferentcharge(TOF_XICs,XICs_R_N,Fract,Fract_cs,TOFAMT_LIST,FractionName,dataname,CS)

if ~isempty(Fract_cs)
    C_ID=TOFAMT_LIST(Fract(Fract_cs),4);
    C_ID_inXICs=Fract(Fract_cs);
    filename=['D:\Program\QTOF_replicate_identification\LongXIC5574f4\Allinteract', num2str(FractionName) ,'_3.ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_0', num2str(dataname),'_f4.centroid1.tolerance20.usemax0.chargeXIC',num2str(CS),'.mat'];
    load(filename);
    eval(['chargeXIC=chargeXIC',num2str(CS),';']);
    for i=1:6
        TOF_XICs{i}(:, C_ID_inXICs)=chargeXIC{i}(1:XICs_R_N, C_ID);
    end
end