function [MZvalue0102_corr,MZvalue0102_non_corr,...
AMZdiff_corr,AMZdiff_non_corr,...
AMASSdiff_corr,AMASSdiff_non_corr,...
AR0102diff_corr,AR0102diff_non_corr,AR0201diff_non_corr,...
AR0103diff_corr,AR0103diff_non_corr,AR0301diff_non_corr,...
AR0203diff_corr,AR0203diff_non_corr,AR0302diff_non_corr,...
AT0102diff_corr,AT0102diff_non_corr,...
AT0201diff_corr,AT0201diff_non_corr,...
AT0103diff_corr,AT0103diff_non_corr,...
AT0301diff_corr,AT0301diff_non_corr,...
AT0203diff_corr,AT0203diff_non_corr,...
AT0302diff_corr,AT0302diff_non_corr,...
PR01diff_corr,PR01diff_non_corr,...
PR02diff_corr,PR02diff_non_corr,...
PR03diff_corr,PR03diff_non_corr]=GenerateAllATARScore(Training_matrix,...
   groundtruthinterval01_final,groundtruthinterval02_final,groundtruthinterval03_final,...
   IntervalList01_final,IntervalList02_final,IntervalList03_final,...
   monoXICs01_final,monoXICs02_final,monoXICs03_final,...
   iso1stXICs01_final,iso1stXICs02_final,iso1stXICs03_final,...
   retentiont01l1,retentiont02l1,retentiont03l1,...
   PP12,PP13,PP21,PP23,PP31,PP32)

MZvalue0102_corr=[];MZvalue0102_non_corr=[];
AMZdiff_corr=[];AMZdiff_non_corr=[];
AMASSdiff_corr=[];AMASSdiff_non_corr=[];

AR0102diff_corr=[];AR0102diff_non_corr=[];AR0201diff_non_corr=[];
AR0103diff_corr=[];AR0103diff_non_corr=[];AR0301diff_non_corr=[];
AR0203diff_corr=[];AR0203diff_non_corr=[];AR0302diff_non_corr=[];
AT0102diff_corr=[];AT0102diff_non_corr=[];
AT0201diff_corr=[];AT0201diff_non_corr=[];
AT0103diff_corr=[];AT0103diff_non_corr=[];
AT0301diff_corr=[];AT0301diff_non_corr=[];
AT0203diff_corr=[];AT0203diff_non_corr=[];
AT0302diff_corr=[];AT0302diff_non_corr=[];

PR01diff_corr=[];PR01diff_non_corr=[];
PR02diff_corr=[];PR02diff_non_corr=[];
PR03diff_corr=[];PR03diff_non_corr=[];


for i=1:size(Training_matrix,1)
    
    id_in01=Training_matrix(i,8);
    id_in02=Training_matrix(i,16);
    id_in03=Training_matrix(i,24);
    
    intervalid01=groundtruthinterval01_final(id_in01,2);
    intervalid02=groundtruthinterval02_final(id_in02,2);
    intervalid03=groundtruthinterval03_final(id_in03,2);

    intervalmatrix01=IntervalList01_final{id_in01}.intervallist;
    intervalmatrix02=IntervalList02_final{id_in02}.intervallist;
    intervalmatrix03=IntervalList03_final{id_in03}.intervallist;
    
    %%%%%% PR 01
    for j=1:size(intervalmatrix01,1)
        scanstart=intervalmatrix01(j,1);
        scanend=intervalmatrix01(j,2);
        monoelutionprofile01=monoXICs01_final(scanstart:scanend,id_in01);
        isoelutionprofile01=iso1stXICs01_final(scanstart:scanend,id_in01);
        newdata_sampling01=monoelutionprofile01';
        newdata_sampling02=isoelutionprofile01';
        [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
        Training_Rstats010203(i).PR01(j)=STATS(1);
    end
    %%%%%%
    %%%%%% PR 02
    for j=1:size(intervalmatrix02,1)
        scanstart=intervalmatrix02(j,1);
        scanend=intervalmatrix02(j,2);
        monoelutionprofile02=monoXICs02_final(scanstart:scanend,id_in02);
        isoelutionprofile02=iso1stXICs02_final(scanstart:scanend,id_in02);
        newdata_sampling01=monoelutionprofile02';
        newdata_sampling02=isoelutionprofile02';
        [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
        Training_Rstats010203(i).PR02(j)=STATS(1);
    end
    %%%%%%
    %%%%%% PR 03
    for j=1:size(intervalmatrix03,1)
        scanstart=intervalmatrix03(j,1);
        scanend=intervalmatrix03(j,2);
        monoelutionprofile03=monoXICs03_final(scanstart:scanend,id_in03);
        isoelutionprofile03=iso1stXICs03_final(scanstart:scanend,id_in03);
        newdata_sampling01=monoelutionprofile03';
        newdata_sampling02=isoelutionprofile03';
        [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
        Training_Rstats010203(i).PR03(j)=STATS(1);
    end
    %%%%%%
    %%%%%% AR 0102 AT 0102 AT0201
    for j1=1:size(intervalmatrix01,1)
        for j2=1:size(intervalmatrix02,1)

                scanstart01=intervalmatrix01(j1,1);
                scanend01=intervalmatrix01(j1,2);
                monoelutionprofile01=monoXICs01_final(scanstart01:scanend01,id_in01);
                isoelutionprofile01=iso1stXICs01_final(scanstart01:scanend01,id_in01);

                scanstart02=intervalmatrix02(j2,1);
                scanend02=intervalmatrix02(j2,2);
                monoelutionprofile02=monoXICs02_final(scanstart02:scanend02,id_in02);
                isoelutionprofile02=iso1stXICs02_final(scanstart02:scanend02,id_in02);

                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01', monoelutionprofile02');
                [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                Training_Rstats010203(i).AR0102(j1,j2)=STATS(1);

                t1=(retentiont01l1(scanstart01)+retentiont01l1(scanend01))/2;
                t2=(retentiont02l1(scanstart02)+retentiont02l1(scanend02))/2;
                t12=polyval(PP12,t1)-t2;
                t21=polyval(PP21,t2)-t1;
                Training_Time010203(i).AT0102(j1,j2)=t12;
                Training_Time010203(i).AT0201(j2,j1)=t21;

        end
    end
    %%%%%% 
    %%%%%% AR 0103 AT 0103 AT0301
    for j1=1:size(intervalmatrix01,1)
        for j3=1:size(intervalmatrix03,1)

                scanstart01=intervalmatrix01(j1,1);
                scanend01=intervalmatrix01(j1,2);
                monoelutionprofile01=monoXICs01_final(scanstart01:scanend01,id_in01);
                isoelutionprofile01=iso1stXICs01_final(scanstart01:scanend01,id_in01);

                scanstart03=intervalmatrix03(j3,1);
                scanend03=intervalmatrix03(j3,2);
                monoelutionprofile03=monoXICs03_final(scanstart03:scanend03,id_in03);
                isoelutionprofile03=iso1stXICs03_final(scanstart03:scanend03,id_in03);

                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01', monoelutionprofile03');
                [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                Training_Rstats010203(i).AR0103(j1,j3)=STATS(1);

                t1=(retentiont01l1(scanstart01)+retentiont01l1(scanend01))/2;
                t3=(retentiont03l1(scanstart03)+retentiont03l1(scanend03))/2;
                t13=polyval(PP13,t1)-t3;
                t31=polyval(PP31,t3)-t1;
                Training_Time010203(i).AT0103(j1,j3)=t13;
                Training_Time010203(i).AT0301(j3,j1)=t31;

        end
    end
    %%%%%%  
    %%%%%% AR 0203 AT 0203 AT0302
    for j2=1:size(intervalmatrix02,1)
        for j3=1:size(intervalmatrix03,1)

                scanstart02=intervalmatrix02(j2,1);
                scanend02=intervalmatrix02(j2,2);
                monoelutionprofile02=monoXICs02_final(scanstart02:scanend02,id_in02);
                isoelutionprofile02=iso1stXICs02_final(scanstart02:scanend02,id_in02);

                scanstart03=intervalmatrix03(j3,1);
                scanend03=intervalmatrix03(j3,2);
                monoelutionprofile03=monoXICs03_final(scanstart03:scanend03,id_in03);
                isoelutionprofile03=iso1stXICs03_final(scanstart03:scanend03,id_in03);

                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile02', monoelutionprofile03');
                [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                Training_Rstats010203(i).AR0203(j2,j3)=STATS(1);

                t2=(retentiont02l1(scanstart02)+retentiont02l1(scanend02))/2;
                t3=(retentiont03l1(scanstart03)+retentiont03l1(scanend03))/2;
                t23=polyval(PP23,t2)-t3;
                t32=polyval(PP32,t3)-t2;
                Training_Time010203(i).AT0203(j2,j3)=t23;
                Training_Time010203(i).AT0302(j3,j2)=t32;

        end
    end
    %%%%%% 
    
    AR0102diff_corr=[AR0102diff_corr,Training_Rstats010203(i).AR0102(intervalid01,intervalid02)];
    AT0102diff_corr=[AT0102diff_corr,Training_Time010203(i).AT0102(intervalid01,intervalid02)];
    AT0201diff_corr=[AT0201diff_corr,Training_Time010203(i).AT0201(intervalid02,intervalid01)];
    
    AR0103diff_corr=[AR0103diff_corr,Training_Rstats010203(i).AR0103(intervalid01,intervalid03)];
    AT0103diff_corr=[AT0103diff_corr,Training_Time010203(i).AT0103(intervalid01,intervalid03)];
    AT0301diff_corr=[AT0301diff_corr,Training_Time010203(i).AT0301(intervalid03,intervalid01)];
    
    AR0203diff_corr=[AR0203diff_corr,Training_Rstats010203(i).AR0203(intervalid02,intervalid03)];
    AT0203diff_corr=[AT0203diff_corr,Training_Time010203(i).AT0203(intervalid02,intervalid03)];
    AT0302diff_corr=[AT0302diff_corr,Training_Time010203(i).AT0302(intervalid03,intervalid02)];
    
      
    PR01diff_corr=[PR01diff_corr,Training_Rstats010203(i).PR01(intervalid01)];
    PR02diff_corr=[PR02diff_corr,Training_Rstats010203(i).PR02(intervalid02)];   
    PR03diff_corr=[PR03diff_corr,Training_Rstats010203(i).PR03(intervalid03)]; 
    
    Prv01=Training_Rstats010203(i).PR01;
    Prv02=Training_Rstats010203(i).PR02;
    Prv03=Training_Rstats010203(i).PR03;
    Prv01(intervalid01)=[];    Prv02(intervalid02)=[];    Prv03(intervalid03)=[];
    PR01diff_non_corr=[PR01diff_non_corr,Prv01];
    PR02diff_non_corr=[PR02diff_non_corr,Prv02];
    PR03diff_non_corr=[PR03diff_non_corr,Prv03];
    
    
    tdiff0102=Training_Time010203(i).AT0102(:,intervalid02);
    tdiff0201=Training_Time010203(i).AT0201(:,intervalid01);
    tdiff0103=Training_Time010203(i).AT0103(:,intervalid03);
    tdiff0301=Training_Time010203(i).AT0301(:,intervalid01);
    tdiff0203=Training_Time010203(i).AT0203(:,intervalid03);
    tdiff0302=Training_Time010203(i).AT0302(:,intervalid02);
    tdiff0102(intervalid01)=[]; tdiff0201(intervalid02)=[];
    tdiff0103(intervalid01)=[]; tdiff0301(intervalid03)=[];
    tdiff0203(intervalid02)=[]; tdiff0302(intervalid03)=[];
    tdiff0102v1=reshape(tdiff0102,1,length(tdiff0102));
    tdiff0201v1=reshape(tdiff0201,1,length(tdiff0201));
    tdiff0103v1=reshape(tdiff0103,1,length(tdiff0103));
    tdiff0301v1=reshape(tdiff0301,1,length(tdiff0301));
    tdiff0203v1=reshape(tdiff0203,1,length(tdiff0203));
    tdiff0302v1=reshape(tdiff0302,1,length(tdiff0302));

    AT0102diff_non_corr=[AT0102diff_non_corr,tdiff0102v1];
    AT0201diff_non_corr=[AT0201diff_non_corr,tdiff0201v1];
    AT0103diff_non_corr=[AT0103diff_non_corr,tdiff0103v1];
    AT0301diff_non_corr=[AT0301diff_non_corr,tdiff0301v1]; 
    AT0203diff_non_corr=[AT0203diff_non_corr,tdiff0203v1]; 
    AT0302diff_non_corr=[AT0302diff_non_corr,tdiff0302v1]; 
    
    
    Rdiff0102=Training_Rstats010203(i).AR0102(:,intervalid02);
    Rdiff0201=Training_Rstats010203(i).AR0102(intervalid01,:);
    Rdiff0103=Training_Rstats010203(i).AR0103(:,intervalid03);
    Rdiff0301=Training_Rstats010203(i).AR0103(intervalid01,:);
    Rdiff0203=Training_Rstats010203(i).AR0203(:,intervalid03);
    Rdiff0302=Training_Rstats010203(i).AR0203(intervalid02,:);
    Rdiff0102(intervalid01)=[]; Rdiff0201(intervalid02)=[];
    Rdiff0103(intervalid01)=[]; Rdiff0301(intervalid03)=[];
    Rdiff0203(intervalid02)=[]; Rdiff0302(intervalid03)=[];
    Rdiff0102v1=reshape(Rdiff0102,1,length(Rdiff0102));
    Rdiff0201v1=reshape(Rdiff0201,1,length(Rdiff0201));
    Rdiff0103v1=reshape(Rdiff0103,1,length(Rdiff0103));
    Rdiff0301v1=reshape(Rdiff0301,1,length(Rdiff0301));
    Rdiff0203v1=reshape(Rdiff0203,1,length(Rdiff0203));
    Rdiff0302v1=reshape(Rdiff0302,1,length(Rdiff0302));
    
    AR0102diff_non_corr=[AR0102diff_non_corr,Rdiff0102v1];
    AR0201diff_non_corr=[AR0201diff_non_corr,Rdiff0201v1];
    AR0103diff_non_corr=[AR0103diff_non_corr,Rdiff0103v1];
    AR0301diff_non_corr=[AR0301diff_non_corr,Rdiff0301v1]; 
    AR0203diff_non_corr=[AR0203diff_non_corr,Rdiff0203v1]; 
    AR0302diff_non_corr=[AR0302diff_non_corr,Rdiff0302v1];    
        
%     Ground_elution_peak01=monoXICs01_final(interval01(1):interval01(2),id_in01);
%     Ground_elution_peak02=monoXICs02_final(interval02(1):interval02(2),id_in02);
%     Ground_elution_peak03=monoXICs03_final(interval03(1):interval03(2),id_in03);   

end

%%%%%%%%%%%%%%%%
