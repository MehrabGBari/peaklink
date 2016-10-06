function data01_XIC=dataXIC_generation_fromTOFQUANTRAW(data01,total_detect_inAMT,ID_5574_01_f4_unique_inAMT,ID_5574_01_f4);

for i=1:length(ID_5574_01_f4_unique_inAMT)
    
    Po_5574_01_f4_unique_inAMT=ID_5574_01_f4_unique_inAMT(i);
    Po_5574_01_f4=find(ID_5574_01_f4==1);
    Po_5574_01_f4_inAMT=find(total_detect_inAMT(Po_5574_01_f4,1)==Po_5574_01_f4_unique_inAMT);
    Po_5574_01_f4_inTOF=total_detect_inAMT(Po_5574_01_f4(Po_5574_01_f4_inAMT),2);  
    
    for j=1:length(Po_5574_01_f4_inTOF)
                
        [peptidenew, massdiffList, isHeavy]=modprocess({data01{Po_5574_01_f4_inTOF(j)}.peptide});
        peptidenew=peptidenew{1};
        [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(data01{Po_5574_01_f4_inTOF(j)}.peptide));
        csvalue=data01{Po_5574_01_f4_inTOF(j)}.charge;
        mzvalue=(weight+csvalue*1.0073)/csvalue;
        total_rt=((data01{Po_5574_01_f4_inTOF(j)}.end_rt+data01{Po_5574_01_f4_inTOF(j)}.start_rt)/2)/data01{Po_5574_01_f4_inTOF(j)}.normal_rt;
        normal_rt=data01{Po_5574_01_f4_inTOF(j)}.normal_rt;
        normal_start_rt=data01{Po_5574_01_f4_inTOF(j)}.start_rt/total_rt;
        normal_end_rt=data01{Po_5574_01_f4_inTOF(j)}.end_rt/total_rt;
        interval=data01{Po_5574_01_f4_inTOF(j)}.interval;
        intervalRawData=data01{Po_5574_01_f4_inTOF(j)}.intervalRawData;
        
%         data01_XIC(i).interval_information(j).mzvalue=mzvalue;
%         data01_XIC(i).interval_information(j).csvalue=csvalue;
%         data01_XIC(i).interval_information(j).total_rt=total_rt;
%         data01_XIC(i).interval_information(j).normal_rt=normal_rt;
%         data01_XIC(i).interval_information(j).normal_start_rt=normal_start_rt;
%         data01_XIC(i).interval_information(j).normal_end_rt=normal_end_rt;
%         data01_XIC(i).interval_information(j).intervalRawData=intervalRawData;
%         data01_XIC(i).interval_information(j).interval=interval;
        
        data01_XIC(i).mzvalue(j)=mzvalue;
        data01_XIC(i).csvalue(j)=csvalue;
        data01_XIC(i).total_rt(j)=total_rt;
        data01_XIC(i).normal_rt(j)=normal_rt;
        data01_XIC(i).normal_start_rt(j)=normal_start_rt;
        data01_XIC(i).normal_end_rt(j)=normal_end_rt;
        data01_XIC(i).intervalRawData{j}=intervalRawData;
        data01_XIC(i).interval{j}=interval;
    
    end
    
end
