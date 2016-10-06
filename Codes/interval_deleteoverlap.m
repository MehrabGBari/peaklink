function data01_ComXICv1=interval_deleteoverlap(data01_ComXIC)

data01_ComXICv1=data01_ComXIC;
for i=1:length(data01_ComXIC)
    Interval=[];
    for m=1:length(data01_ComXIC(i).interval)
        Interval=[Interval;data01_ComXIC(i).interval{m}];
    end
    Intervalv1=Interval;
    num=1;
    while size(Interval,1)>=1
        Scan_start=Interval(1,1);
        Scan_end=Interval(1,2);
        Posi=find((Scan_start-Interval(:,2)).*(Scan_end-Interval(:,1))<=0);
        Raw_ID{num}=Posi;
        Interval_Length=[];
        for k=1:length(Posi)
            Interval_Length=[Interval_Length;Interval(Posi(k),2)-Interval(Posi(k),1)];        
        end
        Raw_Interval_length{num}=Interval_Length;
        num=num+1;
        Interval(Posi,:)=[]; 
    end
    Raw_keep=[];
    intervalRawData=data01_ComXICv1(i).intervalRawData;
    Goodinterval=data01_ComXICv1(i).Goodinterval;
    mzvalue=data01_ComXICv1(i).mzvalue;
    csvalue=data01_ComXICv1(i).csvalue;
    total_rt=data01_ComXICv1(i).total_rt;
    normal_rt=data01_ComXICv1(i).normal_rt;
    normal_start_rt=data01_ComXICv1(i).normal_start_rt;
    normal_end_rt=data01_ComXICv1(i).normal_end_rt;
    interval=data01_ComXICv1(i).interval;
    intervalRawData=data01_ComXICv1(i).intervalRawData;
    iso=data01_ComXICv1(i).iso;
    IndentFileSurName=data01_ComXICv1(i).IndentFileSurName;
    Posi_inTOFQUANTRAW=data01_ComXICv1(i).Posi_inTOFQUANTRAW;    
    
    for j=1:length(Raw_ID)
        [V,P]=max(Raw_Interval_length{j});
        Matrix1{j}=intervalRawData{Raw_ID{j}(P)};
        Matrix2(j)=Goodinterval(Raw_ID{j}(P));
        Matrix3(j)=mzvalue(Raw_ID{j}(P));
        Matrix4(j)=csvalue(Raw_ID{j}(P));
        Matrix5(j)=total_rt(Raw_ID{j}(P));
        Matrix6(j)=normal_rt(Raw_ID{j}(P));
        Matrix7(j)=normal_start_rt(Raw_ID{j}(P));
        Matrix8(j)=normal_end_rt(Raw_ID{j}(P));
        Matrix9{j}=interval{Raw_ID{j}(P)};
        Matrix10{j}=intervalRawData{Raw_ID{j}(P)};
        Matrix11{j}=iso{Raw_ID{j}(P)};
        Matrix12{j}=IndentFileSurName{Raw_ID{j}(P)};
        Matrix13(j)=Posi_inTOFQUANTRAW(Raw_ID{j}(P));
        intervalRawData(Raw_ID{j})=[];
        Goodinterval(Raw_ID{j})=[];
        mzvalue(Raw_ID{j})=[];
        csvalue(Raw_ID{j})=[];
        total_rt(Raw_ID{j})=[];
        normal_rt(Raw_ID{j})=[];
        normal_start_rt(Raw_ID{j})=[];
        normal_end_rt(Raw_ID{j})=[];
        interval(Raw_ID{j})=[];
        iso(Raw_ID{j})=[];
        IndentFileSurName(Raw_ID{j})=[];
        Posi_inTOFQUANTRAW(Raw_ID{j})=[];
    end
    
    
    data01_ComXICv1(i).intervalRawData=Matrix1;
    data01_ComXICv1(i).Goodinterval=Matrix2;
    data01_ComXICv1(i).mzvalue=Matrix3;
    data01_ComXICv1(i).csvalue=Matrix4;
    data01_ComXICv1(i).total_rt=Matrix5;
    data01_ComXICv1(i).normal_rt=Matrix6;
    data01_ComXICv1(i).normal_start_rt=Matrix7;
    data01_ComXICv1(i).normal_end_rt=Matrix8;
    data01_ComXICv1(i).interval=Matrix9;
    data01_ComXICv1(i).intervalRawData=Matrix10;
    data01_ComXICv1(i).iso=Matrix11;
    data01_ComXICv1(i).IndentFileSurName=Matrix12;
    data01_ComXICv1(i).Posi_inTOFQUANTRAW=Matrix13;  
    
    clear Raw_ID Raw_Interval_length Matrix1 Matrix2 Matrix3
    clear Matrix4 Matrix5 Matrix6 Matrix7 Matrix8 Matrix9 Matrix10 Matrix11
    clear Matrix12 Matrix13
end
