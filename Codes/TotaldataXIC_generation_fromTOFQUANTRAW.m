function [data01_TotalXIC_inAMTlist,Empty_ID]=TotaldataXIC_generation_fromTOFQUANTRAW(data01,AMT_database,threshold)

for j=1:length(data01)
    [peptidenew, massdiffList, isHeavy]=modprocess({data01{j}.peptide});
    peptidenew=peptidenew{1};
    [peptideformula,isotopepattern,weight]=aminocalculation_mod(peptidenew,7,getmodificationformula(data01{j}.peptide));
    csvalue=data01{j}.charge;
    mzvalue=(weight+csvalue*1.0073)/csvalue;
    Mzlist_inTOF(j)=mzvalue;
end

Empty_ID=[];
for i=1:length(AMT_database)
    cs=AMT_database(i).charge;
    pepseq=AMT_database(i).peptide;
    mono_mz=AMT_database(i).mono_mz;
    num=1;    
    P_mz=Mzlist_inTOF>=mono_mz-mono_mz*0.00002 & Mzlist_inTOF<=mono_mz+mono_mz*0.00002;
    if sum(P_mz)~=0
        P_mzv1=find(P_mz==1);
        for j=1:length(P_mzv1)
            
            csvalue=data01{P_mzv1(j)}.charge;
            mzvalue=Mzlist_inTOF(P_mzv1(j));
            total_rt=((data01{P_mzv1(j)}.end_rt+data01{P_mzv1(j)}.start_rt)/2)/data01{P_mzv1(j)}.normal_rt;
            normal_rt=data01{P_mzv1(j)}.normal_rt;
            normal_start_rt=data01{P_mzv1(j)}.start_rt/total_rt;
            normal_end_rt=data01{P_mzv1(j)}.end_rt/total_rt;
            interval=data01{P_mzv1(j)}.interval;
            intervalRawData=data01{P_mzv1(j)}.intervalRawData;
            iso=data01{P_mzv1(j)}.iso;
            IndentFileSurName=data01{P_mzv1(j)}.IndentFileSurName;
            
                            %%%%%%%%%%%%%
                            isoList=iso;
                            [~,maxindex]=max(isoList(1:2));
                            temp0=intervalRawData(:,1);
                            temp1=intervalRawData(:,2);
                            temp2=intervalRawData(:,3);
                            temp3=intervalRawData(:,4);
                            temp4=intervalRawData(:,5);
                            temp5=intervalRawData(:,6);
                            temp=[temp0,temp1,temp2,temp3,temp4,temp5];
                            tempmax=temp(:,maxindex);
                            tempmax1=temp(:,3-maxindex);
                            tempmax4=temp(:,maxindex+4);
                            %LC profile should be simliar for first and fifth XIC interval,
                            %first and second; 5th and 6th
                            corr1=cosinecorr(tempmax,tempmax1);
                            corr2=cosinecorr(temp(:,maxindex+4),temp(:,7-maxindex));
                            corr3=cosinecorr(temp(:,maxindex),temp(:,maxindex+4));
                            tcount=0;
                            for nn=1:size(temp,1)
                                if KLdistance(temp0(nn),temp1(nn),isoList(1),isoList(2))<-2.5 && KLdistance(temp2(nn), temp3(nn),sum(temp2), sum(temp3))<-2.5 && KLdistance(temp4(nn),temp5(nn),sum(temp4), sum(temp5))<-2.5 &&  KLdistance(tempmax(nn),tempmax4(nn),sum(tempmax), sum(tempmax4))<-2.5
                                    tcount=tcount+1;
                                end
                            end

                            if corr1^2>threshold && corr2^2>threshold && corr3^2>threshold && KLdistance(isoList(1),isoList(2),sum(temp0), sum(temp1))<-2.5
                                if tcount>4 || 2*tcount>size(temp,1)
                                    data01_TotalXIC_inAMTlist(i).Goodinterval(num)=1;
                                else
                                    data01_TotalXIC_inAMTlist(i).Goodinterval(num)=0;
                                end
                            else
                                data01_TotalXIC_inAMTlist(i).Goodinterval(num)=0;
                            end
                            %%%%%%%%%%%%%
            
            data01_TotalXIC_inAMTlist(i).mzvalue(num)=mzvalue;
            data01_TotalXIC_inAMTlist(i).csvalue(num)=csvalue;
            data01_TotalXIC_inAMTlist(i).total_rt(num)=total_rt;
            data01_TotalXIC_inAMTlist(i).normal_rt(num)=normal_rt;
            data01_TotalXIC_inAMTlist(i).normal_start_rt(num)=normal_start_rt;
            data01_TotalXIC_inAMTlist(i).normal_end_rt(num)=normal_end_rt;
            data01_TotalXIC_inAMTlist(i).interval{num}=interval;
            data01_TotalXIC_inAMTlist(i).intervalRawData{num}=intervalRawData;
            data01_TotalXIC_inAMTlist(i).iso{num}=iso;        
            data01_TotalXIC_inAMTlist(i).IndentFileSurName{num}=IndentFileSurName;
            data01_TotalXIC_inAMTlist(i).Posi_inTOFQUANTRAW(num)=P_mzv1(j);
            num=num+1;
        end
        
    end

    if num==1
        data01_TotalXIC_inAMTlist(i).EmptyXIC=1;
        Empty_ID=[Empty_ID,i];
    end    
end
