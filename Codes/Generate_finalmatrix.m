function Testing_candidate_matrix=Generate_finalmatrix(Testing_candidate_matrix,...
    groundtruthinterval01_final,groundtruthinterval02_final,groundtruthinterval03_final,...
    IntervalList01_final,IntervalList02_final,IntervalList03_final,...
    monoXICs01_final,monoXICs02_final,monoXICs03_final,...
    monoXICs0102_final,monoXICs0103_final,...
    monoXICs0201_final,monoXICs0203_final,...
    monoXICs0301_final,monoXICs0302_final,...
    iso1stXICs01_final,iso1stXICs02_final,iso1stXICs03_final,...
    retentiont01l1,retentiont02l1,retentiont03l1,...
    AR0102_PARMHAT1,AR0102_PARMHAT2,AR0201_PARMHAT1,AR0201_PARMHAT2,...
    AR0103_PARMHAT1,AR0103_PARMHAT2,AR0301_PARMHAT1,AR0301_PARMHAT2,...
    AR0203_PARMHAT1,AR0203_PARMHAT2,AR0302_PARMHAT1,AR0302_PARMHAT2,...
    timemodel0102_MU_sig,timemodel0102_SIGMA_sig,...
    timemodel0201_MU_sig,timemodel0201_SIGMA_sig,...
    timemodel0103_MU_sig,timemodel0103_SIGMA_sig,...
    timemodel0301_MU_sig,timemodel0301_SIGMA_sig,...
    timemodel0203_MU_sig,timemodel0203_SIGMA_sig,...
    timemodel0302_MU_sig,timemodel0302_SIGMA_sig,...
    timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig,...
    timemodel0201_MU_notsig,timemodel0201_SIGMA_notsig,...
    timemodel0103_MU_notsig,timemodel0103_SIGMA_notsig,...
    timemodel0301_MU_notsig,timemodel0301_SIGMA_notsig,...
    timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig,...
    timemodel0302_MU_notsig,timemodel0302_SIGMA_notsig,...
    PP12,PP13,PP21,PP23,PP31,PP32)



Record_testingmatrix=zeros(size(Testing_candidate_matrix,1),7);
for i=1:size(Testing_candidate_matrix,1)
    
    Index01=Testing_candidate_matrix(i,8);
    Index02=Testing_candidate_matrix(i,16);
    Index03=Testing_candidate_matrix(i,24);
    
    Prob01=Testing_candidate_matrix(i,7);
    Prob02=Testing_candidate_matrix(i,15);
    Prob03=Testing_candidate_matrix(i,23);
    Record_testingmatrix(i,4:6)=[Prob01,Prob02,Prob03];
    
    Jud_vector=[Index01,Index02,Index03]~=0;
    Record_testingmatrix(i,1:3)=Jud_vector;
    Jud_num=Jud_vector*[1 2 4]';
    switch Jud_num
        
        case 1  %%%% [1 0 0]         
                XIC01=monoXICs01_final(:,Index01);
                intervalid01=groundtruthinterval01_final(Index01,2);
                intervalmatrix01=IntervalList01_final{Index01}.intervallist;
                Ground_interval01=intervalmatrix01(intervalid01,:);
                ElutionPeak01=XIC01(Ground_interval01(1):Ground_interval01(2));
                
                XIC0102=monoXICs0102_final(:,Index01);
                XIC0103=monoXICs0103_final(:,Index01);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0102]=intervaldetection(XIC0102,mininterval,maxnointervals);
                [intervalmatrix0103]=intervaldetection(XIC0103,mininterval,maxnointervals);
                
                for j=1:size(intervalmatrix0102,1)                    
                    Scanstart0102=intervalmatrix0102(j,1);
                    Scanend0102=intervalmatrix0102(j,2);
                    if Scanend0102~=0 && Scanstart0102~=0
                            ElutionPeak0102=XIC0102(Scanstart0102:Scanend0102);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak01', ElutionPeak0102');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0102(j)=STATS(1);

                            Testing_Rstats010203(i).AR0102_corrProb(j)=gampdf(1-STATS(1),AR0102_PARMHAT1(1),AR0102_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0102_noncorrProb(j)=gampdf(STATS(1),AR0102_PARMHAT2(1),AR0102_PARMHAT2(2));                   

                            t1=(retentiont01l1(Ground_interval01(1))+retentiont01l1(Ground_interval01(2)))/2;
                            t2=(retentiont02l1(Scanstart0102)+retentiont02l1(Scanend0102))/2;
                            t=polyval(PP12,t1)-t2;
                            Testing_Time010203(i).AT0102(j)=t;
                            Testing_Time010203(i).AT0102_corrProb(j)=normpdf(t,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
                            Testing_Time010203(i).AT0102_noncorrProb(j)=normpdf(t,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);                   
                    else
                            Testing_Rstats010203(i).AR0102(j)=0;
                            Testing_Rstats010203(i).AR0102_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0102_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0102(j)=1000000;
                            Testing_Time010203(i).AT0102_corrProb(j)=0;
                            Testing_Time010203(i).AT0102_noncorrProb(j)=0.1;                  
        
                    end
                end
                [V,P]=max(log(Testing_Rstats010203(i).AR0102_corrProb)+log(Testing_Time010203(i).AT0102_corrProb));
                Testing_candidate_matrix(i,13:14)=intervalmatrix0102(P,:);
                
                for j=1:size(intervalmatrix0103,1)                    
                    Scanstart0103=intervalmatrix0103(j,1);
                    Scanend0103=intervalmatrix0103(j,2);
                    if Scanstart0103~=0 && Scanend0103~=0 && Scanend0103-Scanstart0103<=7500
                            ElutionPeak0103=XIC0103(Scanstart0103:Scanend0103);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak01', ElutionPeak0103');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0103(j)=STATS(1);

                            Testing_Rstats010203(i).AR0103_corrProb(j)=gampdf(1-STATS(1),AR0103_PARMHAT1(1),AR0103_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0103_noncorrProb(j)=gampdf(STATS(1),AR0103_PARMHAT2(1),AR0103_PARMHAT2(2));                   

                            t1=(retentiont01l1(Ground_interval01(1))+retentiont01l1(Ground_interval01(2)))/2;
                            t2=(retentiont03l1(Scanstart0103)+retentiont03l1(Scanend0103))/2;
                            t=polyval(PP13,t1)-t2;
                            Testing_Time010203(i).AT0103(j)=t;
                            Testing_Time010203(i).AT0103_corrProb(j)=normpdf(t,timemodel0103_MU_sig,timemodel0103_SIGMA_sig);
                            Testing_Time010203(i).AT0103_noncorrProb(j)=normpdf(t,timemodel0103_MU_notsig,timemodel0103_SIGMA_notsig);                   
                    else 
                            Testing_Rstats010203(i).AR0103(j)=0;
                            Testing_Rstats010203(i).AR0103_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0103_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0103(j)=1000000;
                            Testing_Time010203(i).AT0103_corrProb(j)=0;
                            Testing_Time010203(i).AT0103_noncorrProb(j)=0.1;  
                        
                    end
                end
                [V,P]=max(log(Testing_Rstats010203(i).AR0103_corrProb)+log(Testing_Time010203(i).AT0103_corrProb));
                Testing_candidate_matrix(i,21:22)=intervalmatrix0103(P,:);
                
        case 2  %%%% [0 1 0]

                XIC02=monoXICs02_final(:,Index02);
                intervalid02=groundtruthinterval02_final(Index02,2);
                intervalmatrix02=IntervalList02_final{Index02}.intervallist;
                Ground_interval02=intervalmatrix02(intervalid02,:);
                ElutionPeak02=XIC02(Ground_interval02(1):Ground_interval02(2));
                
                XIC0201=monoXICs0201_final(:,Index02);
                XIC0203=monoXICs0203_final(:,Index02);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0201]=intervaldetection(XIC0201,mininterval,maxnointervals);
                [intervalmatrix0203]=intervaldetection(XIC0203,mininterval,maxnointervals);
                
                for j=1:size(intervalmatrix0201,1)                    
                    Scanstart0201=intervalmatrix0201(j,1);
                    Scanend0201=intervalmatrix0201(j,2);
                    if Scanend0201~=0 && Scanstart0201~=0
                        if Scanend0201-Scanstart0201<=7500
                            ElutionPeak0201=XIC0201(Scanstart0201:Scanend0201);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak02', ElutionPeak0201');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0201(j)=STATS(1);
                        else Testing_Rstats010203(i).AR0201(j)=0.1;
                        end

                            Testing_Rstats010203(i).AR0201_corrProb(j)=gampdf(1-STATS(1),AR0201_PARMHAT1(1),AR0201_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0201_noncorrProb(j)=gampdf(STATS(1),AR0201_PARMHAT2(1),AR0201_PARMHAT2(2));                   

                            t1=(retentiont02l1(Ground_interval02(1))+retentiont02l1(Ground_interval02(2)))/2;
                            t2=(retentiont01l1(Scanstart0201)+retentiont01l1(Scanend0201))/2;
                            t=polyval(PP21,t1)-t2;
                            Testing_Time010203(i).AT0201(j)=t;
                            Testing_Time010203(i).AT0201_corrProb(j)=normpdf(t,timemodel0201_MU_sig,timemodel0201_SIGMA_sig);
                            Testing_Time010203(i).AT0201_noncorrProb(j)=normpdf(t,timemodel0201_MU_notsig,timemodel0201_SIGMA_notsig);                   
                    else
                            Testing_Rstats010203(i).AR0201(j)=0;
                            Testing_Rstats010203(i).AR0201_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0201_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0201(j)=1000000;
                            Testing_Time010203(i).AT0201_corrProb(j)=0;
                            Testing_Time010203(i).AT0201_noncorrProb(j)=0.1;                  
        
                    end
                end
                [V,P]=max(log(Testing_Rstats010203(i).AR0201_corrProb)+log(Testing_Time010203(i).AT0201_corrProb));
                Testing_candidate_matrix(i,5:6)=intervalmatrix0201(P,:);
                
                for j=1:size(intervalmatrix0203,1)                    
                    Scanstart0203=intervalmatrix0203(j,1);
                    Scanend0203=intervalmatrix0203(j,2);
                    if Scanstart0203~=0 && Scanend0203~=0
                            ElutionPeak0203=XIC0203(Scanstart0203:Scanend0203);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak02', ElutionPeak0203');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0203(j)=STATS(1);

                            Testing_Rstats010203(i).AR0203_corrProb(j)=gampdf(1-STATS(1),AR0203_PARMHAT1(1),AR0203_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0203_noncorrProb(j)=gampdf(STATS(1),AR0203_PARMHAT2(1),AR0203_PARMHAT2(2));                   

                            t1=(retentiont02l1(Ground_interval02(1))+retentiont02l1(Ground_interval02(2)))/2;
                            t2=(retentiont03l1(Scanstart0203)+retentiont03l1(Scanend0203))/2;
                            t=polyval(PP23,t1)-t2;
                            Testing_Time010203(i).AT0203(j)=t;
                            Testing_Time010203(i).AT0203_corrProb(j)=normpdf(t,timemodel0203_MU_sig,timemodel0203_SIGMA_sig);
                            Testing_Time010203(i).AT0203_noncorrProb(j)=normpdf(t,timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig);                   
                    else 
                            Testing_Rstats010203(i).AR0203(j)=0;
                            Testing_Rstats010203(i).AR0203_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0203_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0203(j)=1000000;
                            Testing_Time010203(i).AT0203_corrProb(j)=0;
                            Testing_Time010203(i).AT0203_noncorrProb(j)=0.1;  
                        
                    end
                end
                [V,P]=max(log(Testing_Rstats010203(i).AR0203_corrProb)+log(Testing_Time010203(i).AT0203_corrProb));
                Testing_candidate_matrix(i,21:22)=intervalmatrix0203(P,:);
            
        
        case 4  %%%% [0 0 1]
        
            
                XIC03=monoXICs03_final(:,Index03);
                intervalid03=groundtruthinterval03_final(Index03,2);
                intervalmatrix03=IntervalList03_final{Index03}.intervallist;
                Ground_interval03=intervalmatrix03(intervalid03,:);
                ElutionPeak03=XIC03(Ground_interval03(1):Ground_interval03(2));
                
                XIC0301=monoXICs0301_final(:,Index03);
                XIC0302=monoXICs0302_final(:,Index03);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0301]=intervaldetection(XIC0301,mininterval,maxnointervals);
                [intervalmatrix0302]=intervaldetection(XIC0302,mininterval,maxnointervals);
                
                for j=1:size(intervalmatrix0301,1)                    
                    Scanstart0301=intervalmatrix0301(j,1);
                    Scanend0301=intervalmatrix0301(j,2);
                    if Scanend0301~=0 && Scanstart0301~=0
                            ElutionPeak0301=XIC0301(Scanstart0301:Scanend0301);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak03', ElutionPeak0301');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0301(j)=STATS(1);

                            Testing_Rstats010203(i).AR0301_corrProb(j)=gampdf(1-STATS(1),AR0301_PARMHAT1(1),AR0301_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0301_noncorrProb(j)=gampdf(STATS(1),AR0301_PARMHAT2(1),AR0301_PARMHAT2(2));                   

                            t1=(retentiont03l1(Ground_interval03(1))+retentiont03l1(Ground_interval03(2)))/2;
                            t2=(retentiont01l1(Scanstart0301)+retentiont01l1(Scanend0301))/2;
                            t=polyval(PP31,t1)-t2;
                            Testing_Time010203(i).AT0301(j)=t;
                            Testing_Time010203(i).AT0301_corrProb(j)=normpdf(t,timemodel0301_MU_sig,timemodel0301_SIGMA_sig);
                            Testing_Time010203(i).AT0301_noncorrProb(j)=normpdf(t,timemodel0301_MU_notsig,timemodel0301_SIGMA_notsig);                   
                    else
                            Testing_Rstats010203(i).AR0301(j)=0;
                            Testing_Rstats010203(i).AR0301_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0301_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0301(j)=1000000;
                            Testing_Time010203(i).AT0301_corrProb(j)=0;
                            Testing_Time010203(i).AT0301_noncorrProb(j)=0.1;                  
        
                    end
                end
                [V,P]=max(log(Testing_Rstats010203(i).AR0301_corrProb)+log(Testing_Time010203(i).AT0301_corrProb));
                Testing_candidate_matrix(i,5:6)=intervalmatrix0301(P,:);
                
                for j=1:size(intervalmatrix0302,1)                    
                    Scanstart0302=intervalmatrix0302(j,1);
                    Scanend0302=intervalmatrix0302(j,2);
                    if Scanstart0302~=0 && Scanend0302~=0
                            ElutionPeak0302=XIC0302(Scanstart0302:Scanend0302);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak03', ElutionPeak0302');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0302(j)=STATS(1);

                            Testing_Rstats010203(i).AR0302_corrProb(j)=gampdf(1-STATS(1),AR0302_PARMHAT1(1),AR0302_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0302_noncorrProb(j)=gampdf(STATS(1),AR0302_PARMHAT2(1),AR0302_PARMHAT2(2));                   

                            t1=(retentiont03l1(Ground_interval03(1))+retentiont03l1(Ground_interval03(2)))/2;
                            t2=(retentiont02l1(Scanstart0302)+retentiont02l1(Scanend0302))/2;
                            t=polyval(PP32,t1)-t2;
                            Testing_Time010203(i).AT0302(j)=t;
                            Testing_Time010203(i).AT0302_corrProb(j)=normpdf(t,timemodel0302_MU_sig,timemodel0302_SIGMA_sig);
                            Testing_Time010203(i).AT0302_noncorrProb(j)=normpdf(t,timemodel0302_MU_notsig,timemodel0302_SIGMA_notsig);                   
                    else 
                            Testing_Rstats010203(i).AR0302(j)=0;
                            Testing_Rstats010203(i).AR0302_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0302_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0302(j)=1000000;
                            Testing_Time010203(i).AT0302_corrProb(j)=0;
                            Testing_Time010203(i).AT0302_noncorrProb(j)=0.1;  
                        
                    end
                end
                [V,P]=max(log(Testing_Rstats010203(i).AR0302_corrProb)+log(Testing_Time010203(i).AT0302_corrProb));
                Testing_candidate_matrix(i,13:14)=intervalmatrix0302(P,:);
            
            
            
        case 3  %%%% [1 1 0]

                XIC01=monoXICs01_final(:,Index01);
                intervalid01=groundtruthinterval01_final(Index01,2);
                intervalmatrix01=IntervalList01_final{Index01}.intervallist;
                Ground_interval01=intervalmatrix01(intervalid01,:);
                ElutionPeak01=XIC01(Ground_interval01(1):Ground_interval01(2));
                
                XIC0103=monoXICs0103_final(:,Index01);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0103]=intervaldetection(XIC0103,mininterval,maxnointervals);
            
                XIC02=monoXICs02_final(:,Index02);
                intervalid02=groundtruthinterval02_final(Index02,2);
                intervalmatrix02=IntervalList02_final{Index02}.intervallist;
                Ground_interval02=intervalmatrix02(intervalid02,:);
                ElutionPeak02=XIC02(Ground_interval02(1):Ground_interval02(2));
                
                XIC0203=monoXICs0203_final(:,Index02);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0203]=intervaldetection(XIC0203,mininterval,maxnointervals);
            
                for j=1:size(intervalmatrix0103,1)                    
                    Scanstart0103=intervalmatrix0103(j,1);
                    Scanend0103=intervalmatrix0103(j,2);
                    if Scanstart0103~=0 && Scanend0103~=0
                            ElutionPeak0103=XIC0103(Scanstart0103:Scanend0103);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak01', ElutionPeak0103');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0103(j)=STATS(1);

                            Testing_Rstats010203(i).AR0103_corrProb(j)=gampdf(1-STATS(1),AR0103_PARMHAT1(1),AR0103_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0103_noncorrProb(j)=gampdf(STATS(1),AR0103_PARMHAT2(1),AR0103_PARMHAT2(2));                   

                            t1=(retentiont01l1(Ground_interval01(1))+retentiont01l1(Ground_interval01(2)))/2;
                            t2=(retentiont03l1(Scanstart0103)+retentiont03l1(Scanend0103))/2;
                            t=polyval(PP13,t1)-t2;
                            Testing_Time010203(i).AT0103(j)=t;
                            Testing_Time010203(i).AT0103_corrProb(j)=normpdf(t,timemodel0103_MU_sig,timemodel0103_SIGMA_sig);
                            Testing_Time010203(i).AT0103_noncorrProb(j)=normpdf(t,timemodel0103_MU_notsig,timemodel0103_SIGMA_notsig);                   
                    else 
                            Testing_Rstats010203(i).AR0103(j)=0;
                            Testing_Rstats010203(i).AR0103_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0103_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0103(j)=1000000;
                            Testing_Time010203(i).AT0103_corrProb(j)=0;
                            Testing_Time010203(i).AT0103_noncorrProb(j)=0.1;  
                        
                    end
                end
                [V0103,P0103]=max(log(Testing_Rstats010203(i).AR0103_corrProb)+log(Testing_Time010203(i).AT0103_corrProb));
            
                for j=1:size(intervalmatrix0203,1)                    
                    Scanstart0203=intervalmatrix0203(j,1);
                    Scanend0203=intervalmatrix0203(j,2);
                    if Scanstart0203~=0 && Scanend0203~=0
                            ElutionPeak0203=XIC0203(Scanstart0203:Scanend0203);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak02', ElutionPeak0203');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0203(j)=STATS(1);

                            Testing_Rstats010203(i).AR0203_corrProb(j)=gampdf(1-STATS(1),AR0203_PARMHAT1(1),AR0203_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0203_noncorrProb(j)=gampdf(STATS(1),AR0203_PARMHAT2(1),AR0203_PARMHAT2(2));                   

                            t1=(retentiont02l1(Ground_interval02(1))+retentiont02l1(Ground_interval02(2)))/2;
                            t2=(retentiont03l1(Scanstart0203)+retentiont03l1(Scanend0203))/2;
                            t=polyval(PP23,t1)-t2;
                            Testing_Time010203(i).AT0203(j)=t;
                            Testing_Time010203(i).AT0203_corrProb(j)=normpdf(t,timemodel0203_MU_sig,timemodel0203_SIGMA_sig);
                            Testing_Time010203(i).AT0203_noncorrProb(j)=normpdf(t,timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig);                   
                    else 
                            Testing_Rstats010203(i).AR0203(j)=0;
                            Testing_Rstats010203(i).AR0203_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0203_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0203(j)=1000000;
                            Testing_Time010203(i).AT0203_corrProb(j)=0;
                            Testing_Time010203(i).AT0203_noncorrProb(j)=0.1;  
                        
                    end
                end
                [V0203,P0203]=max(log(Testing_Rstats010203(i).AR0203_corrProb)+log(Testing_Time010203(i).AT0203_corrProb));
                
                if V0103>=V0203
                    Testing_candidate_matrix(i,21:22)=intervalmatrix0103(P0103,:);
                else 
                    Testing_candidate_matrix(i,21:22)=intervalmatrix0203(P0203,:);
                end
                
                if P0103==P0203
                     Record_testingmatrix(i,7)=1;
                end
         
        case 5  %%%% [1 0 1]

                XIC01=monoXICs01_final(:,Index01);
                intervalid01=groundtruthinterval01_final(Index01,2);
                intervalmatrix01=IntervalList01_final{Index01}.intervallist;
                Ground_interval01=intervalmatrix01(intervalid01,:);
                ElutionPeak01=XIC01(Ground_interval01(1):Ground_interval01(2));
                
                XIC0102=monoXICs0102_final(:,Index01);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0102]=intervaldetection(XIC0102,mininterval,maxnointervals);
         
                XIC03=monoXICs03_final(:,Index03);
                intervalid03=groundtruthinterval03_final(Index03,2);
                intervalmatrix03=IntervalList03_final{Index03}.intervallist;
                Ground_interval03=intervalmatrix03(intervalid03,:);
                ElutionPeak03=XIC03(Ground_interval03(1):Ground_interval03(2));
                
                XIC0302=monoXICs0302_final(:,Index03);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0302]=intervaldetection(XIC0302,mininterval,maxnointervals);
                
                for j=1:size(intervalmatrix0102,1)                    
                    Scanstart0102=intervalmatrix0102(j,1);
                    Scanend0102=intervalmatrix0102(j,2);
                    if Scanend0102~=0 && Scanstart0102~=0
                            ElutionPeak0102=XIC0102(Scanstart0102:Scanend0102);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak01', ElutionPeak0102');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0102(j)=STATS(1);

                            Testing_Rstats010203(i).AR0102_corrProb(j)=gampdf(1-STATS(1),AR0102_PARMHAT1(1),AR0102_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0102_noncorrProb(j)=gampdf(STATS(1),AR0102_PARMHAT2(1),AR0102_PARMHAT2(2));                   

                            t1=(retentiont01l1(Ground_interval01(1))+retentiont01l1(Ground_interval01(2)))/2;
                            t2=(retentiont02l1(Scanstart0102)+retentiont02l1(Scanend0102))/2;
                            t=polyval(PP12,t1)-t2;
                            Testing_Time010203(i).AT0102(j)=t;
                            Testing_Time010203(i).AT0102_corrProb(j)=normpdf(t,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
                            Testing_Time010203(i).AT0102_noncorrProb(j)=normpdf(t,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);                   
                    else
                            Testing_Rstats010203(i).AR0102(j)=0;
                            Testing_Rstats010203(i).AR0102_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0102_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0102(j)=1000000;
                            Testing_Time010203(i).AT0102_corrProb(j)=0;
                            Testing_Time010203(i).AT0102_noncorrProb(j)=0.1;                  
        
                    end
                end
                [V0102,P0102]=max(log(Testing_Rstats010203(i).AR0102_corrProb)+log(Testing_Time010203(i).AT0102_corrProb));
                
                for j=1:size(intervalmatrix0302,1)                    
                    Scanstart0302=intervalmatrix0302(j,1);
                    Scanend0302=intervalmatrix0302(j,2);
                    if Scanstart0302~=0 && Scanend0302~=0
                            ElutionPeak0302=XIC0302(Scanstart0302:Scanend0302);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak03', ElutionPeak0302');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0302(j)=STATS(1);

                            Testing_Rstats010203(i).AR0302_corrProb(j)=gampdf(1-STATS(1),AR0302_PARMHAT1(1),AR0302_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0302_noncorrProb(j)=gampdf(STATS(1),AR0302_PARMHAT2(1),AR0302_PARMHAT2(2));                   

                            t1=(retentiont03l1(Ground_interval03(1))+retentiont03l1(Ground_interval03(2)))/2;
                            t2=(retentiont02l1(Scanstart0302)+retentiont02l1(Scanend0302))/2;
                            t=polyval(PP32,t1)-t2;
                            Testing_Time010203(i).AT0302(j)=t;
                            Testing_Time010203(i).AT0302_corrProb(j)=normpdf(t,timemodel0302_MU_sig,timemodel0302_SIGMA_sig);
                            Testing_Time010203(i).AT0302_noncorrProb(j)=normpdf(t,timemodel0302_MU_notsig,timemodel0302_SIGMA_notsig);                   
                    else 
                            Testing_Rstats010203(i).AR0302(j)=0;
                            Testing_Rstats010203(i).AR0302_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0302_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0302(j)=1000000;
                            Testing_Time010203(i).AT0302_corrProb(j)=0;
                            Testing_Time010203(i).AT0302_noncorrProb(j)=0.1;  
                        
                    end
                end
                [V0302,P0302]=max(log(Testing_Rstats010203(i).AR0302_corrProb)+log(Testing_Time010203(i).AT0302_corrProb));
                
                if V0102>=V0302
                        Testing_candidate_matrix(i,13:14)=intervalmatrix0102(P0102,:);
                else
                        Testing_candidate_matrix(i,13:14)=intervalmatrix0302(P0302,:);
                end
                
                if P0102==P0302
                     Record_testingmatrix(i,7)=1;
                end
            
        case 6  %%%% [0 1 1]
            
                XIC02=monoXICs02_final(:,Index02);
                intervalid02=groundtruthinterval02_final(Index02,2);
                intervalmatrix02=IntervalList02_final{Index02}.intervallist;
                Ground_interval02=intervalmatrix02(intervalid02,:);
                ElutionPeak02=XIC02(Ground_interval02(1):Ground_interval02(2));
                
                XIC0201=monoXICs0201_final(:,Index02);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0201]=intervaldetection(XIC0201,mininterval,maxnointervals);
            
                XIC03=monoXICs03_final(:,Index03);
                intervalid03=groundtruthinterval03_final(Index03,2);
                intervalmatrix03=IntervalList03_final{Index03}.intervallist;
                Ground_interval03=intervalmatrix03(intervalid03,:);
                ElutionPeak03=XIC03(Ground_interval03(1):Ground_interval03(2));
                
                XIC0301=monoXICs0301_final(:,Index03);
                mininterval=4; maxnointervals=6;
                [intervalmatrix0301]=intervaldetection(XIC0301,mininterval,maxnointervals);

                for j=1:size(intervalmatrix0201,1)                    
                    Scanstart0201=intervalmatrix0201(j,1);
                    Scanend0201=intervalmatrix0201(j,2);
                    if Scanend0201-Scanstart0201>=8000
                        Scanend0201=0; Scanstart0201=0;
                    end
                    if Scanend0201~=0 && Scanstart0201~=0
                            ElutionPeak0201=XIC0201(Scanstart0201:Scanend0201);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak02', ElutionPeak0201');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0201(j)=STATS(1);

                            Testing_Rstats010203(i).AR0201_corrProb(j)=gampdf(1-STATS(1),AR0201_PARMHAT1(1),AR0201_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0201_noncorrProb(j)=gampdf(STATS(1),AR0201_PARMHAT2(1),AR0201_PARMHAT2(2));                   

                            t1=(retentiont02l1(Ground_interval02(1))+retentiont02l1(Ground_interval02(2)))/2;
                            t2=(retentiont01l1(Scanstart0201)+retentiont01l1(Scanend0201))/2;
                            t=polyval(PP21,t1)-t2;
                            Testing_Time010203(i).AT0201(j)=t;
                            Testing_Time010203(i).AT0201_corrProb(j)=normpdf(t,timemodel0201_MU_sig,timemodel0201_SIGMA_sig);
                            Testing_Time010203(i).AT0201_noncorrProb(j)=normpdf(t,timemodel0201_MU_notsig,timemodel0201_SIGMA_notsig);                   

                    else
                            Testing_Rstats010203(i).AR0201(j)=0;
                            Testing_Rstats010203(i).AR0201_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0201_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0201(j)=1000000;
                            Testing_Time010203(i).AT0201_corrProb(j)=0;
                            Testing_Time010203(i).AT0201_noncorrProb(j)=0.1;                  
        
                    end
                end
                [V0201,P0201]=max(log(Testing_Rstats010203(i).AR0201_corrProb)+log(Testing_Time010203(i).AT0201_corrProb));


                for j=1:size(intervalmatrix0301,1)                    
                    Scanstart0301=intervalmatrix0301(j,1);
                    Scanend0301=intervalmatrix0301(j,2);
                    if Scanend0301-Scanstart0301>=8000
                        Scanend0301=0;Scanstart0301=0;
                    end
                    if Scanend0301~=0 && Scanstart0301~=0
                            ElutionPeak0301=XIC0301(Scanstart0301:Scanend0301);
                            [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak03', ElutionPeak0301');
                            [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                            Testing_Rstats010203(i).AR0301(j)=STATS(1);

                            Testing_Rstats010203(i).AR0301_corrProb(j)=gampdf(1-STATS(1),AR0301_PARMHAT1(1),AR0301_PARMHAT1(2));
                            Testing_Rstats010203(i).AR0301_noncorrProb(j)=gampdf(STATS(1),AR0301_PARMHAT2(1),AR0301_PARMHAT2(2));                   

                            t1=(retentiont03l1(Ground_interval03(1))+retentiont03l1(Ground_interval03(2)))/2;
                            t2=(retentiont01l1(Scanstart0301)+retentiont01l1(Scanend0301))/2;
                            t=polyval(PP31,t1)-t2;
                            Testing_Time010203(i).AT0301(j)=t;
                            Testing_Time010203(i).AT0301_corrProb(j)=normpdf(t,timemodel0301_MU_sig,timemodel0301_SIGMA_sig);
                            Testing_Time010203(i).AT0301_noncorrProb(j)=normpdf(t,timemodel0301_MU_notsig,timemodel0301_SIGMA_notsig);                   
                    else
                            Testing_Rstats010203(i).AR0301(j)=0;
                            Testing_Rstats010203(i).AR0301_corrProb(j)=0;
                            Testing_Rstats010203(i).AR0301_noncorrProb(j)=0.1;                   

                            Testing_Time010203(i).AT0301(j)=1000000;
                            Testing_Time010203(i).AT0301_corrProb(j)=0;
                            Testing_Time010203(i).AT0301_noncorrProb(j)=0.1;                  
        
                    end
                end
                [V0301,P0301]=max(log(Testing_Rstats010203(i).AR0301_corrProb)+log(Testing_Time010203(i).AT0301_corrProb));
               
                if V0201>=V0301
                        Testing_candidate_matrix(i,5:6)=intervalmatrix0201(P0201,:);
                else                
                        Testing_candidate_matrix(i,5:6)=intervalmatrix0301(P0301,:);
                end
                
                if P0201==P0301
                     Record_testingmatrix(i,7)=1;
                end
             
    end

end