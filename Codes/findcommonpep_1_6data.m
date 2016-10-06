clc
clear all

load data01_information
load data02_information
load data03_information
% load data04_information
% load data05_information
% load data06_information

data01_information.remarks=1*ones(length(data01_information.peptide),1);
data02_information.remarks=2*ones(length(data02_information.peptide),1);
data03_information.remarks=3*ones(length(data03_information.peptide),1);
% data04_information.remarks=4*ones(length(data04_information.peptide),1);
% data05_information.remarks=5*ones(length(data05_information.peptide),1);
% data06_information.remarks=6*ones(length(data06_information.peptide),1);

data01v1=filteroverlapp(data01_information);
data02v1=filteroverlapp(data02_information);
data03v1=filteroverlapp(data03_information);
% data04v1=filteroverlapp(data04_information);
% data05v1=filteroverlapp(data05_information);
% data06v1=filteroverlapp(data06_information);

P_threshold=0.95;
Id01_95=find(data01v1.probability>=P_threshold);
Id02_95=find(data02v1.probability>=P_threshold);
Id03_95=find(data03v1.probability>=P_threshold);
% Id04_95=find(data04v1.probability>=P_threshold);
% Id05_95=find(data05v1.probability>=P_threshold);
% Id06_95=find(data06v1.probability>=P_threshold);


% fprintf('\n Please choose the input information xls file: \n')
% pepinf=xls2struct();

pep01_original=data01v1.peptide;
chargestate01_original=data01v1.assumed_charge;
ms2information01_original=data01v1.retention_time_sec;
mass01_original=data01v1.calc_neutral_pep_mass;
remarks01_original=data01v1.remarks;
Probability01_original=data01v1.probability;
mzratio01_original=data01v1.MZratio;
mzdiff01_original=(mass01_original./chargestate01_original+1.0073)-mzratio01_original;



pep02_original=data02v1.peptide;
chargestate02_original=data02v1.assumed_charge;
ms2information02_original=data02v1.retention_time_sec;
mass02_original=data02v1.calc_neutral_pep_mass;
remarks02_original=data02v1.remarks;
Probability02_original=data02v1.probability;
mzratio02_original=data02v1.MZratio;

pep03_original=data03v1.peptide;
chargestate03_original=data03v1.assumed_charge;
ms2information03_original=data03v1.retention_time_sec;
mass03_original=data03v1.calc_neutral_pep_mass;
remarks03_original=data03v1.remarks;
Probability03_original=data03v1.probability;
mzratio03_original=data03v1.MZratio;
% clear pepinf

% [retentiont01l1,data01l1raw,data01l1,retentiont01_l1l2,MZInt01_l1l2]=readrawdata();
% [retentiont02l1,data02l1raw,data02l1,retentiont02_l1l2,MZInt02_l1l2]=readrawdata();
% [retentiont03l1,data03l1raw,data03l1,retentiont03_l1l2,MZInt03_l1l2]=readrawdata();

% save data01levelone data01l1  retentiont01l1 retentiont01_l1l2 MZInt01_l1l2
% save data02levelone data02l1  retentiont02l1 retentiont02_l1l2 MZInt02_l1l2
% save data03levelone data03l1  retentiont03l1 retentiont03_l1l2 MZInt03_l1l2
load data01levelone
load data02levelone
load data03levelone
%%%%%%% the retention time format is different so all are changed to be one
%%%%%%% row
% retentiont01l1;%%% row
% retentiont02l1=retentiont02l1(:,1)';%%% row
% retentiont03l1=retentiont03l1';%%% row


%%%%%%%%%%%%%%%%%%%%%%%
pep01=pep01_original(Id01_95);
chargestate01=chargestate01_original(Id01_95);
ms2information01=ms2information01_original(Id01_95);
mass01=mass01_original(Id01_95);
remarks01=remarks01_original(Id01_95);
Probability01=Probability01_original(Id01_95);

pep02=pep02_original(Id02_95);
chargestate02=chargestate02_original(Id02_95);
ms2information02=ms2information02_original(Id02_95);
mass02=mass02_original(Id02_95);
remarks02=remarks02_original(Id02_95);
Probability02=Probability02_original(Id02_95);

pep03=pep03_original(Id03_95);
chargestate03=chargestate03_original(Id03_95);
ms2information03=ms2information03_original(Id03_95);
mass03=mass03_original(Id03_95);
remarks03=remarks03_original(Id03_95);
Probability03=Probability03_original(Id03_95);


% for i=1:length(pep01)
%     peplist01_finalv1{i}=pep01{i}(3:end-2);
% end
% for i=1:length(pep02)
%     peplist02_finalv1{i}=pep02{i}(3:end-2);
% end
% 
% Num1=0;
% for i=1:length(pep01)
%             TF=strcmp(peplist01_finalv1{i},peplist02_finalv1);
%             Posi=find(TF==1);
%             if isempty(Posi)==0
%                 Num1=Num1+1;
%             end
% end

% %%%%%%%% get XICs
% [monoXICs01,iso1stXICs01]=getpepXICmatrix(data01l1, pep01, chargestate01,mass01);
% save data01_XICs monoXICs01 iso1stXICs01 -v7.3
% [monoXICs0102,iso1stXICs0102]=getpepXICmatrix(data02l1, pep01, chargestate01,mass01);
% save data0102_XICs monoXICs0102 iso1stXICs0102 -v7.3
% [monoXICs0103,iso1stXICs0103]=getpepXICmatrix(data03l1, pep01, chargestate01,mass01);
% save data0103_XICs monoXICs0103 iso1stXICs0103 -v7.3
% 
% [monoXICs02,iso1stXICs02]=getpepXICmatrix(data02l1, pep02, chargestate02,mass02);
% save data02_XICs monoXICs02 iso1stXICs02 -v7.3
% [monoXICs0201,iso1stXICs0201]=getpepXICmatrix(data01l1, pep02, chargestate02,mass02);
% save data0201_XICs monoXICs0201 iso1stXICs0201 -v7.3
% [monoXICs0203,iso1stXICs0203]=getpepXICmatrix(data03l1, pep02, chargestate02,mass02);
% save data0203_XICs monoXICs0203 iso1stXICs0203 -v7.3
% 
% [monoXICs03,iso1stXICs03]=getpepXICmatrix(data03l1, pep03, chargestate03,mass03);
% save data03_XICs monoXICs03 iso1stXICs03 -v7.3
% [monoXICs0301,iso1stXICs0301]=getpepXICmatrix(data01l1, pep03, chargestate03,mass03);
% save data0301_XICs monoXICs0301 iso1stXICs0301 -v7.3
% [monoXICs0302,iso1stXICs0302]=getpepXICmatrix(data02l1, pep03, chargestate03,mass03);
% save data0302_XICs monoXICs0302 iso1stXICs0302 -v7.3

load data01_XICs
load data0102_XICs
load data0103_XICs
load data02_XICs
load data0201_XICs
load data0203_XICs
load data03_XICs
load data0301_XICs
load data0302_XICs


% %%%%%%%%% verify pep
% [peplist01_final, chargestate01_final, groundtruthinterval01_final, IntervalList01_final, ms2time01_final, monoXICs01_final, iso1stXICs01_final]=verifypeptide(pep01,chargestate01,ms2information01,monoXICs01,iso1stXICs01,retentiont01l1);
% [peplist02_final, chargestate02_final, groundtruthinterval02_final, IntervalList02_final, ms2time02_final, monoXICs02_final, iso1stXICs02_final]=verifypeptide(pep02,chargestate02,ms2information02,monoXICs02,iso1stXICs02,retentiont02l1);
% [peplist03_final, chargestate03_final, groundtruthinterval03_final, IntervalList03_final, ms2time03_final, monoXICs03_final, iso1stXICs03_final]=verifypeptide(pep03,chargestate03,ms2information03,monoXICs03,iso1stXICs03,retentiont03l1);
% % %%%%%%%%%
% save data01demoinformation monoXICs01_final iso1stXICs01_final peplist01_final chargestate01_final groundtruthinterval01_final IntervalList01_final ms2time01_final
% save data02demoinformation monoXICs02_final iso1stXICs02_final peplist02_final chargestate02_final groundtruthinterval02_final IntervalList02_final ms2time02_final
% save data03demoinformation monoXICs03_final iso1stXICs03_final peplist03_final chargestate03_final groundtruthinterval03_final IntervalList03_final ms2time03_final

load data01demoinformation
load data02demoinformation
load data03demoinformation

monoXICs0102_final=monoXICs0102(:,groundtruthinterval01_final(:,1));
monoXICs0103_final=monoXICs0103(:,groundtruthinterval01_final(:,1));
monoXICs0201_final=monoXICs0201(:,groundtruthinterval02_final(:,1));
monoXICs0203_final=monoXICs0203(:,groundtruthinterval02_final(:,1));
monoXICs0301_final=monoXICs0301(:,groundtruthinterval03_final(:,1));
monoXICs0302_final=monoXICs0302(:,groundtruthinterval03_final(:,1));

%%%%%%%%%%
Interval01_ground=[];
Interval02_ground=[];
Interval03_ground=[];
for i=1:length(IntervalList01_final)
    Interval01_ground=[Interval01_ground;IntervalList01_final{i}.intervallist(groundtruthinterval01_final(i,2),:)];
end
for i=1:length(IntervalList02_final)
    Interval02_ground=[Interval02_ground;IntervalList02_final{i}.intervallist(groundtruthinterval02_final(i,2),:)];
end
for i=1:length(IntervalList03_final)
    Interval03_ground=[Interval03_ground;IntervalList03_final{i}.intervallist(groundtruthinterval03_final(i,2),:)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Total.peptide=[peplist01_final;peplist02_final;peplist03_final];
Total.assumed_charge=[chargestate01_final;chargestate02_final;chargestate03_final];
Total.timepoint=[ms2time01_final';ms2time02_final';ms2time03_final'];
posi01=groundtruthinterval01_final(:,1);
posi02=groundtruthinterval02_final(:,1);
posi03=groundtruthinterval03_final(:,1);
Total.mass=[mass01(posi01);mass02(posi02);mass03(posi03)];
Total.remarks=[remarks01(posi01);remarks02(posi02);remarks03(posi03)];
Total.probability=[Probability01(posi01);Probability02(posi02);Probability03(posi03)];
Total.interval=[Interval01_ground;Interval02_ground;Interval03_ground];
Total.Indexineachdata=[[1:length(peplist01_final)]';[1:length(peplist02_final)]';[1:length(peplist03_final)]'];
Total.beta=zeros(length(Total.peptide),1);

save Total Total

Num=1;
Totalinformationmatrix=[];
for i=1:length(Total.peptide)
    if Total.beta(i)==0
            Pepseq{Num}=Total.peptide{i};
            Num=Num+1;
            TF=strcmp(Total.peptide{i},Total.peptide);
            Posi=find(TF==1);
            inf_vector=zeros(1,24);
            if length(Posi)>1
                for k=1:length(Posi);
                    Index_data=Total.remarks(Posi(k));
                    inf_vector((Index_data-1)*8+1:Index_data*8)=[Total.assumed_charge(Posi(k)), Total.mass(Posi(k)), Total.timepoint(Posi(k)), Total.remarks(Posi(k)), Total.interval(Posi(k),:), Total.probability(Posi(k)), Total.Indexineachdata(Posi(k))];
                end        
                Totalinformationmatrix=[Totalinformationmatrix;inf_vector];
                Total.beta(Posi)=1;
            else
                Index_data=Total.remarks(Posi);
                inf_vector((Index_data-1)*8+1:Index_data*8)=[Total.assumed_charge(Posi), Total.mass(Posi), Total.timepoint(Posi), Total.remarks(Posi), Total.interval(Posi,:), Total.probability(Posi), Total.Indexineachdata(Posi)];
                Totalinformationmatrix=[Totalinformationmatrix;inf_vector];    
                Total.beta(Posi)=1;
            end    
    end
end


%%%%%%%%%%%% The column of Totalinformationmatrix are:
%%%%%%%%%%%% cs01 mass01 timepoint01 remarks01(from data X) intervalstart01
%%%%%%%%%%%% intervalend01 Prob01 indexindata01 (1:8) 
%%%%%%%%%%%% cs02 mass02 timepoint02 remarks02(from data X) intervalstart02
%%%%%%%%%%%% intervalend02 Prob02 indexindata02 (9:16) 
%%%%%%%%%%%% cs03 mass03 timepoint03 remarks03(from data X) intervalstart03
%%%%%%%%%%%% intervalend03 Prob03 indexindata03 (17:24)
%%%%%%%%%%%%
Remarks_Matrix=[Totalinformationmatrix(:,4),Totalinformationmatrix(:,12),Totalinformationmatrix(:,20)];
sum(Remarks_Matrix,1)./[1 2 3];

Remark_PlusTotal=sum(Remarks_Matrix,2);
Remark_Plus12=Remarks_Matrix(:,1)+Remarks_Matrix(:,2);
Remark_Plus13=Remarks_Matrix(:,1)+Remarks_Matrix(:,3);
Remark_Plus23=Remarks_Matrix(:,2)+Remarks_Matrix(:,3);

IDOverlap_total=find(Remark_PlusTotal==6);
IDOverlap_12=find(Remark_Plus12==3);
IDOverlap_13=find(Remark_Plus13==4);
IDOverlap_23=find(Remark_Plus23==5);

Training_candidate_matrix=Totalinformationmatrix(IDOverlap_total,:);
Testing_candidate_matrix_original=Totalinformationmatrix;
Testing_candidate_matrix_original(IDOverlap_total,:)=[];
Testing_Pepseq_original=Pepseq;
Testing_Pepseq_original(IDOverlap_total)=[];

Id_CSnotsame_intesting=[];
for i=1:size(Testing_candidate_matrix_original,1)
    
    CS_vector=Testing_candidate_matrix_original(i,[1,9,17]);
    ID_csnull=CS_vector==0;
    Jud_num=ID_csnull*[1,2,4]';
    switch Jud_num
        
        case 1 %(0 1 1) for cs valuess and (1 0 0) for zeros
            if CS_vector(2)~=CS_vector(3)
                Id_CSnotsame_intesting=[Id_CSnotsame_intesting; i];
            end            
        case 2 %(1 0 1) for cs valuess and (0 1 0) for zeros
            if CS_vector(1)~=CS_vector(3)
                Id_CSnotsame_intesting=[Id_CSnotsame_intesting; i];                
            end                
        case 4 %(1 1 0) for cs valuess and (0 0 1) for zeros
            if CS_vector(1)~=CS_vector(2)
                Id_CSnotsame_intesting=[Id_CSnotsame_intesting; i];
            end                
    end
    
end

Testing_candidate_matrix=Testing_candidate_matrix_original;
Testing_candidate_matrix(Id_CSnotsame_intesting,:)=[];
Testing_Pepseq=Testing_Pepseq_original;
Testing_Pepseq(Id_CSnotsame_intesting)=[];
%%%%%%%%%%%%%%% find same charge state pep in 01 02 03
ID_samecs_trainingcandidate=[];
for i=1:size(Training_candidate_matrix,1)
      if Training_candidate_matrix(i,1)==Training_candidate_matrix(i,9) && Training_candidate_matrix(i,1)==Training_candidate_matrix(i,17)
           ID_samecs_trainingcandidate=[ID_samecs_trainingcandidate;i];         
      end
end
Training_candidate_matrixv1=Training_candidate_matrix(ID_samecs_trainingcandidate,:);
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% find pep with ground truth peaks high than 10^6
Training_matrix=[];
for i=1:size(Training_candidate_matrixv1,1)
    
    id_in01=Training_candidate_matrixv1(i,8);
    id_in02=Training_candidate_matrixv1(i,16);
    id_in03=Training_candidate_matrixv1(i,24);
    
    id_intervalindex01=groundtruthinterval01_final(id_in01,2);
    id_intervalindex02=groundtruthinterval02_final(id_in02,2);
    id_intervalindex03=groundtruthinterval03_final(id_in03,2);
    
    interval01=IntervalList01_final{id_in01}.intervallist(id_intervalindex01,:);
    interval02=IntervalList02_final{id_in02}.intervallist(id_intervalindex02,:);
    interval03=IntervalList03_final{id_in03}.intervallist(id_intervalindex03,:);
    
    Ground_elution_peak01=monoXICs01_final(interval01(1):interval01(2),id_in01);
    Ground_elution_peak02=monoXICs02_final(interval02(1):interval02(2),id_in02);
    Ground_elution_peak03=monoXICs03_final(interval03(1):interval03(2),id_in03);
    
    highest_peak01=max(Ground_elution_peak01);
    highest_peak02=max(Ground_elution_peak02);
    highest_peak03=max(Ground_elution_peak03);
    
    if  highest_peak01>=1000000 && highest_peak02>=1000000 && highest_peak03>=1000000
        Training_matrix=[Training_matrix;Training_candidate_matrixv1(i,:)];        
    end
    
    
end
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% training data AT AR value for model

for i=1:size(Training_matrix,1)
    
    id_in01=Training_matrix(i,8);
    id_in02=Training_matrix(i,16);
    id_in03=Training_matrix(i,24);
    
    id_intervalindex01=groundtruthinterval01_final(id_in01,2);
    id_intervalindex02=groundtruthinterval02_final(id_in02,2);
    id_intervalindex03=groundtruthinterval03_final(id_in03,2);
    
    interval01=IntervalList01_final{id_in01}.intervallist(id_intervalindex01,:);
    interval02=IntervalList02_final{id_in02}.intervallist(id_intervalindex02,:);
    interval03=IntervalList03_final{id_in03}.intervallist(id_intervalindex03,:);
    
    Ground_elution_peak01=monoXICs01_final(interval01(1):interval01(2),id_in01);
    Ground_elution_peak02=monoXICs02_final(interval02(1):interval02(2),id_in02);
    Ground_elution_peak03=monoXICs03_final(interval03(1):interval03(2),id_in03);
    
    %%%%%% retentiont format is different
    Ground_trainingtimepoint01(i)=sum(retentiont01l1(interval01))/2;
    Ground_trainingtimepoint02(i)=sum(retentiont02l1(interval02))/2;
    Ground_trainingtimepoint03(i)=sum(retentiont03l1(interval03))/2;
    
    
end

PP12=polyfit(Ground_trainingtimepoint01,Ground_trainingtimepoint02,4);
PP13=polyfit(Ground_trainingtimepoint01,Ground_trainingtimepoint03,4);
PP21=polyfit(Ground_trainingtimepoint02,Ground_trainingtimepoint01,4);
PP23=polyfit(Ground_trainingtimepoint02,Ground_trainingtimepoint03,4);
PP31=polyfit(Ground_trainingtimepoint03,Ground_trainingtimepoint01,4);
PP32=polyfit(Ground_trainingtimepoint03,Ground_trainingtimepoint02,4);
% figure
% plot(Ground_trainingtimepoint01,Ground_trainingtimepoint03,'r.')
% XXX=1:50:14000; YYY=polyval(PP13,XXX);
% hold on; plot(XXX,YYY);


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

%%%%%%%%%%%%%%% get parameter
for i=1:1
    [timemodel0102_MU_sig,timemodel0102_SIGMA_sig] = normfit(AT0102diff_corr);
    [timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig] = normfit(AT0102diff_non_corr);
    [timemodel0201_MU_sig,timemodel0201_SIGMA_sig] = normfit(AT0201diff_corr);
    [timemodel0201_MU_notsig,timemodel0201_SIGMA_notsig] = normfit(AT0201diff_non_corr);
    [timemodel0103_MU_sig,timemodel0103_SIGMA_sig] = normfit(AT0103diff_corr);
    [timemodel0103_MU_notsig,timemodel0103_SIGMA_notsig] = normfit(AT0103diff_non_corr);
    [timemodel0301_MU_sig,timemodel0301_SIGMA_sig] = normfit(AT0301diff_corr);
    [timemodel0301_MU_notsig,timemodel0301_SIGMA_notsig] = normfit(AT0301diff_non_corr);
    [timemodel0203_MU_sig,timemodel0203_SIGMA_sig] = normfit(AT0203diff_corr);
    [timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig] = normfit(AT0203diff_non_corr);
    [timemodel0302_MU_sig,timemodel0302_SIGMA_sig] = normfit(AT0302diff_corr);
    [timemodel0302_MU_notsig,timemodel0302_SIGMA_notsig] = normfit(AT0302diff_non_corr);

    aa= PR01diff_non_corr>=0 & PR01diff_non_corr<=1;
    PR01diff_non_corr_noNaN=PR01diff_non_corr(aa);
    aa= PR01diff_corr>=0 & PR01diff_corr<=1;
    PR01diff_corr_noNaN=PR01diff_corr(aa);
    PR01_PARMHAT1=gamfit(PR01diff_corr_noNaN);
    PR01_PARMHAT2=gamfit(PR01diff_non_corr_noNaN);

    aa= PR02diff_non_corr>=0 & PR02diff_non_corr<=1;
    PR02diff_non_corr_noNaN=PR02diff_non_corr(aa);
    aa= PR02diff_corr>=0 & PR02diff_corr<=1;
    PR02diff_corr_noNaN=PR02diff_corr(aa);
    PR02_PARMHAT1=gamfit(PR02diff_corr_noNaN);
    PR02_PARMHAT2=gamfit(PR02diff_non_corr_noNaN);
    
    aa= PR03diff_non_corr>=0 & PR03diff_non_corr<=1;
    PR03diff_non_corr_noNaN=PR03diff_non_corr(aa);
    aa= PR03diff_corr>=0 & PR03diff_corr<=1;
    PR03diff_corr_noNaN=PR03diff_corr(aa);
    PR03_PARMHAT1=gamfit(PR03diff_corr_noNaN);
    PR03_PARMHAT2=gamfit(PR03diff_non_corr_noNaN);

    aa= AR0102diff_corr>=0 & AR0102diff_corr<=1;
    AR0102diff_corr_noNaN=AR0102diff_corr(aa);
    aa= AR0102diff_non_corr>=0 & AR0102diff_non_corr<=1;
    AR0102diff_non_corr_noNaN=AR0102diff_non_corr(aa);
    aa= AR0201diff_non_corr>=0 & AR0201diff_non_corr<=1;
    AR0201diff_non_corr_noNaN=AR0201diff_non_corr(aa);
    AR0102_PARMHAT1=gamfit(1-AR0102diff_corr_noNaN);
    AR0102_PARMHAT2=gamfit(AR0102diff_non_corr_noNaN);
    AR0201_PARMHAT1=AR0102_PARMHAT1;
    AR0201_PARMHAT2=gamfit(AR0201diff_non_corr_noNaN);
    
    aa= AR0103diff_corr>=0 & AR0103diff_corr<=1;
    AR0103diff_corr_noNaN=AR0103diff_corr(aa);
    aa= AR0103diff_non_corr>=0 & AR0103diff_non_corr<=1;
    AR0103diff_non_corr_noNaN=AR0103diff_non_corr(aa);
    aa= AR0301diff_non_corr>=0 & AR0301diff_non_corr<=1;
    AR0301diff_non_corr_noNaN=AR0301diff_non_corr(aa);
    AR0103_PARMHAT1=gamfit(1-AR0103diff_corr_noNaN);
    AR0103_PARMHAT2=gamfit(AR0103diff_non_corr_noNaN);
    AR0301_PARMHAT1=AR0103_PARMHAT1;
    AR0301_PARMHAT2=gamfit(AR0301diff_non_corr_noNaN);
    
    aa= AR0203diff_corr>=0 & AR0203diff_corr<=1;
    AR0203diff_corr_noNaN=AR0203diff_corr(aa);
    aa= AR0203diff_non_corr>=0 & AR0203diff_non_corr<=1;
    AR0203diff_non_corr_noNaN=AR0203diff_non_corr(aa);
    aa= AR0302diff_non_corr>=0 & AR0302diff_non_corr<=1;
    AR0302diff_non_corr_noNaN=AR0302diff_non_corr(aa);
    AR0203_PARMHAT1=gamfit(1-AR0203diff_corr_noNaN);
    AR0203_PARMHAT2=gamfit(AR0203diff_non_corr_noNaN);
    AR0302_PARMHAT1=AR0203_PARMHAT1;
    AR0302_PARMHAT2=gamfit(AR0302diff_non_corr_noNaN);
    
end
save Parameters AR0102_PARMHAT1 AR0102_PARMHAT2 AR0201_PARMHAT1 AR0201_PARMHAT2 ...
                            AR0103_PARMHAT1 AR0103_PARMHAT2 AR0301_PARMHAT1 AR0301_PARMHAT2 ...
                            AR0203_PARMHAT1 AR0203_PARMHAT2 AR0302_PARMHAT1 AR0302_PARMHAT2 ...
                            timemodel0102_MU_sig timemodel0102_SIGMA_sig ...
                            timemodel0201_MU_sig timemodel0201_SIGMA_sig ...
                            timemodel0103_MU_sig timemodel0103_SIGMA_sig ...
                            timemodel0301_MU_sig timemodel0301_SIGMA_sig ...
                            timemodel0203_MU_sig timemodel0203_SIGMA_sig ...
                            timemodel0302_MU_sig timemodel0302_SIGMA_sig ...
                            timemodel0102_MU_notsig timemodel0102_SIGMA_notsig ...
                            timemodel0201_MU_notsig timemodel0201_SIGMA_notsig ...
                            timemodel0103_MU_notsig timemodel0103_SIGMA_notsig ...
                            timemodel0301_MU_notsig timemodel0301_SIGMA_notsig ...
                            timemodel0203_MU_notsig timemodel0203_SIGMA_notsig ...
                            timemodel0302_MU_notsig timemodel0302_SIGMA_notsig ...
                            PP12 PP21 PP13 PP31 PP23 PP32

%%%%%%%%%%%%%%%

% %%%%%%%%%%% pic
% for i=1:1
%             XX=0:0.1:1;
%             [a,b]=hist(AR0102diff_corr,XX);
%             [c,d]=hist(AR0102diff_non_corr,XX);
%             [e,f]=hist(AR0201diff_non_corr,XX);
%             b2=0:0.05:1;
%             Y0102_sig= gampdf(1-b2,AR0102_PARMHAT1(1),AR0102_PARMHAT1(2));
%             Y0102_notsig= gampdf(b2,AR0102_PARMHAT2(1),AR0102_PARMHAT2(2));
%             Y0201_notsig= gampdf(b2,AR0201_PARMHAT2(1),AR0201_PARMHAT2(2));
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;
%             plot(d,c/sum(c));plot(f,e/sum(e),'k');
%             legend('corr AR0102 ','non corr AR0102','non corr AR0201')
%             subplot(2,1,2);plot(b2,Y0102_sig,'r');hold on;
%             plot(b2,Y0102_notsig);plot(b2,Y0201_notsig,'k')
% 
%             XX=0:0.1:1;
%             [a,b]=hist(AR0103diff_corr,XX);
%             [c,d]=hist(AR0103diff_non_corr,XX);
%             [e,f]=hist(AR0301diff_non_corr,XX);
%             b2=0:0.05:1;
%             Y0103_sig= gampdf(1-b2,AR0103_PARMHAT1(1),AR0103_PARMHAT1(2));
%             Y0103_notsig= gampdf(b2,AR0103_PARMHAT2(1),AR0103_PARMHAT2(2));
%             Y0301_notsig= gampdf(b2,AR0301_PARMHAT2(1),AR0301_PARMHAT2(2));
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;
%             plot(d,c/sum(c));plot(f,e/sum(e),'k');
%             legend('corr AR0103 ','non corr AR0103','non corr AR0301')
%             subplot(2,1,2);plot(b2,Y0103_sig,'r');hold on;
%             plot(b2,Y0103_notsig);plot(b2,Y0301_notsig,'k')
% 
%             XX=0:0.1:1;
%             [a,b]=hist(AR0203diff_corr,XX);
%             [c,d]=hist(AR0203diff_non_corr,XX);
%             [e,f]=hist(AR0302diff_non_corr,XX);
%             b2=0:0.05:1;
%             Y0203_sig= gampdf(1-b2,AR0203_PARMHAT1(1),AR0203_PARMHAT1(2));
%             Y0203_notsig= gampdf(b2,AR0203_PARMHAT2(1),AR0203_PARMHAT2(2));
%             Y0302_notsig= gampdf(b2,AR0302_PARMHAT2(1),AR0302_PARMHAT2(2));
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;
%             plot(d,c/sum(c));plot(f,e/sum(e),'k');
%             legend('corr AR0203 ','non corr AR0203','non corr AR0302')
%             subplot(2,1,2);plot(b2,Y0203_sig,'r');hold on;
%             plot(b2,Y0203_notsig);plot(b2,Y0302_notsig,'k')
% 
%             XX=-12000:100:12000;
%             [a,b]=hist(AT0102diff_corr,XX);
%             [c,d]=hist(AT0102diff_non_corr,XX);
%             t=-12000:100:12000;
%             P_timenorm_sigv2=normpdf(t,timemodel0102_MU_sig,timemodel0102_SIGMA_sig);
%             P_timenorm_notsigv2=normpdf(t,timemodel0102_MU_notsig,timemodel0102_SIGMA_notsig);
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;plot(d,c/sum(c))
%             legend('corr AT0102','non corr AT0102')
%             subplot(2,1,2);plot(t,P_timenorm_sigv2,'r');hold on;plot(t,P_timenorm_notsigv2);
% 
%             XX=-12000:100:12000;
%             [a,b]=hist(AT0201diff_corr,XX);
%             [c,d]=hist(AT0201diff_non_corr,XX);
%             t=-12000:100:12000;
%             P_timenorm_sigv2=normpdf(t,timemodel0201_MU_sig,timemodel0201_SIGMA_sig);
%             P_timenorm_notsigv2=normpdf(t,timemodel0201_MU_notsig,timemodel0201_SIGMA_notsig);
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;plot(d,c/sum(c))
%             legend('corr AT0201','non corr AT0201')
%             subplot(2,1,2);plot(t,P_timenorm_sigv2,'r');hold on;plot(t,P_timenorm_notsigv2);
% 
%             XX=-12000:100:12000;
%             [a,b]=hist(AT0103diff_corr,XX);
%             [c,d]=hist(AT0103diff_non_corr,XX);
%             t=-12000:100:12000;
%             P_timenorm_sigv2=normpdf(t,timemodel0103_MU_sig,timemodel0103_SIGMA_sig);
%             P_timenorm_notsigv2=normpdf(t,timemodel0103_MU_notsig,timemodel0103_SIGMA_notsig);
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;plot(d,c/sum(c))
%             legend('corr AT0103','non corr AT0103')
%             subplot(2,1,2);plot(t,P_timenorm_sigv2,'r');hold on;plot(t,P_timenorm_notsigv2);
% 
%             XX=-12000:100:12000;
%             [a,b]=hist(AT0301diff_corr,XX);
%             [c,d]=hist(AT0301diff_non_corr,XX);
%             t=-12000:100:12000;
%             P_timenorm_sigv2=normpdf(t,timemodel0301_MU_sig,timemodel0301_SIGMA_sig);
%             P_timenorm_notsigv2=normpdf(t,timemodel0301_MU_notsig,timemodel0301_SIGMA_notsig);
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;plot(d,c/sum(c))
%             legend('corr AT0301','non corr AT0301')
%             subplot(2,1,2);plot(t,P_timenorm_sigv2,'r');hold on;plot(t,P_timenorm_notsigv2);
% 
%             XX=-12000:100:12000;
%             [a,b]=hist(AT0203diff_corr,XX);
%             [c,d]=hist(AT0203diff_non_corr,XX);
%             t=-12000:100:12000;
%             P_timenorm_sigv2=normpdf(t,timemodel0203_MU_sig,timemodel0203_SIGMA_sig);
%             P_timenorm_notsigv2=normpdf(t,timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig);
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;plot(d,c/sum(c))
%             legend('corr AT0203','non corr AT0203')
%             subplot(2,1,2);plot(t,P_timenorm_sigv2,'r');hold on;plot(t,P_timenorm_notsigv2);
% 
%             XX=-12000:100:12000;
%             [a,b]=hist(AT0302diff_corr,XX);
%             [c,d]=hist(AT0302diff_non_corr,XX);
%             t=-12000:100:12000;
%             P_timenorm_sigv2=normpdf(t,timemodel0302_MU_sig,timemodel0302_SIGMA_sig);
%             P_timenorm_notsigv2=normpdf(t,timemodel0302_MU_notsig,timemodel0302_SIGMA_notsig);
%             figure
%             subplot(2,1,1);plot(b,a/sum(a),'r');hold on;plot(d,c/sum(c))
%             legend('corr AT0302','non corr AT0302')
%             subplot(2,1,2);plot(t,P_timenorm_sigv2,'r');hold on;plot(t,P_timenorm_notsigv2);
% end
% %%%%%%%%%%%

%%%%%%%%%%%% Get the score for other peptides(not training)
%%%%%%%%%%%% The column of Testing_candidate_matrix are:
%%%%%%%%%%%% cs01 mass01 timepoint01 remarks01(from data X) intervalstart01
%%%%%%%%%%%% intervalend01 Prob01 indexindata01 (1:8)
%%%%%%%%%%%% cs02 mass02 timepoint02 remarks02(from data X) intervalstart02
%%%%%%%%%%%% intervalend02 Prob02 indexindata02 (9:16)
%%%%%%%%%%%% cs02 mass03 timepoint03 remarks03(from data X) intervalstart03
%%%%%%%%%%%% intervalend03 Prob03 indexindata03 (17:24)

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

for i=1:size(Record_testingmatrix,1)
    Prob_max_vector(i)=max(Record_testingmatrix(i,4:6));
end
% Id_Prob_99=find(Prob_max_vector>=0.99);
Id_Prob_99=Prob_max_vector>=0.99;
sum(Id_Prob_99)

id_app2twice=sum(Record_testingmatrix(:,1:3),2)==2;
id_app2once=sum(Record_testingmatrix(:,1:3),2)==1;
Jud_same=Record_testingmatrix(id_app2twice,7);
Id_Jud_same=find(Jud_same==1);

ID_Data01_detect=Testing_candidate_matrix(:,6)~=0;
ID_Data02_detect=Testing_candidate_matrix(:,14)~=0;
ID_Data03_detect=Testing_candidate_matrix(:,22)~=0;
Jud_matrix=[ID_Data01_detect ID_Data02_detect ID_Data03_detect];
ID_3=find(sum(Jud_matrix,2)==3);
ID_2=find(sum(Jud_matrix,2)==2);
ID_1=find(sum(Jud_matrix,2)==1);

ID_3_T=sum(Jud_matrix,2)==3;
ID_3_T_Prob99=Id_Prob_99'+ID_3_T==2;
sum(ID_3_T)/length(ID_3_T)
sum(ID_3_T_Prob99)/sum(Id_Prob_99)
%%%%%%%%%%%%

matrix{1,1}='pepseq';
matrix{1,2}='cs01';matrix{1,3}='mass01';matrix{1,4}='timepoint01';matrix{1,5}='remarks01';
matrix{1,6}='INTstart01';matrix{1,7}='INTend01';matrix{1,8}='Prob01';matrix{1,9}='INDindata01';
matrix{1,10}='cs02';matrix{1,11}='mass02';matrix{1,12}='timepoint02';matrix{1,13}='remarks02';
matrix{1,14}='INTstart02';matrix{1,15}='INTend02';matrix{1,16}='Prob02';matrix{1,17}='INDindata02';
matrix{1,18}='cs03';matrix{1,19}='mass03';matrix{1,20}='timepoint03';matrix{1,21}='remarks03';
matrix{1,22}='INTstart03';matrix{1,23}='INTend03';matrix{1,24}='Prob03';matrix{1,25}='INDindata03';

for i=1:size(Testing_candidate_matrix,1)
    matrix{i+1,1}=Testing_Pepseq{i};
    matrix{i+1,2}=num2str(Testing_candidate_matrix(i,1));
    matrix{i+1,3}=num2str(Testing_candidate_matrix(i,2));
    matrix{i+1,4}=num2str(Testing_candidate_matrix(i,3));
    matrix{i+1,5}=num2str(Testing_candidate_matrix(i,4));
    matrix{i+1,6}=num2str(Testing_candidate_matrix(i,5));
    matrix{i+1,7}=num2str(Testing_candidate_matrix(i,6));
    matrix{i+1,8}=num2str(Testing_candidate_matrix(i,7));
    matrix{i+1,9}=num2str(Testing_candidate_matrix(i,8));
    matrix{i+1,10}=num2str(Testing_candidate_matrix(i,9));
    matrix{i+1,11}=num2str(Testing_candidate_matrix(i,10));
    matrix{i+1,12}=num2str(Testing_candidate_matrix(i,11));
    matrix{i+1,13}=num2str(Testing_candidate_matrix(i,12));
    matrix{i+1,14}=num2str(Testing_candidate_matrix(i,13));
    matrix{i+1,15}=num2str(Testing_candidate_matrix(i,14));
    matrix{i+1,16}=num2str(Testing_candidate_matrix(i,15));
    matrix{i+1,17}=num2str(Testing_candidate_matrix(i,16));
    matrix{i+1,18}=num2str(Testing_candidate_matrix(i,17));
    matrix{i+1,19}=num2str(Testing_candidate_matrix(i,18));
    matrix{i+1,20}=num2str(Testing_candidate_matrix(i,19));
    matrix{i+1,21}=num2str(Testing_candidate_matrix(i,20));
    matrix{i+1,22}=num2str(Testing_candidate_matrix(i,21));
    matrix{i+1,23}=num2str(Testing_candidate_matrix(i,22));
    matrix{i+1,24}=num2str(Testing_candidate_matrix(i,23));
    matrix{i+1,25}=num2str(Testing_candidate_matrix(i,24));
end

FILE=[input('\n Write the output file name: \n', 's')];
xlswrite(FILE,matrix);

% FILE=['..\output_excel\',input('\n Write the output file name: \n', 's')];
% xlswrite(FILE,matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%% verify the SCIFA algorithm
ID_TF02=Totalinformationmatrix(:,16)~=0;
ID_TF03=Totalinformationmatrix(:,24)~=0;
ID_TF0203_verify=find(ID_TF02+ID_TF03==2);
length(ID_TF0203_verify)

Total_verify_matrix=Totalinformationmatrix(ID_TF0203_verify,:);

%%%%%%%%%%%% cs02 mass02 timepoint01 remarks01(from data X) intervalstart01
%%%%%%%%%%%% intervalend01 Prob01 indexindata01 (1:8)
%%%%%%%%%%%% cs02 mass02 timepoint02 remarks02(from data X) intervalstart02
%%%%%%%%%%%% intervalend02 Prob02 indexindata02 (9:16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% peplist01_final      3891                              peplist02_final  5280
%%% chargestate01_final                                  chargestate02_final
%%% groundtruthinterval01_final                     groundtruthinterval02_final
%%% IntervalList01_final                                   IntervalList02_final
%%% ms2time01_final                                         ms2time02_final
%%% monoXICs01_final                                     monoXICs02_final
%%% iso1stXICs01_final                                     iso1stXICs02_final
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% find same charge state pep in 01 02 03
ID_samecs_0203verify=[];
for i=1:size(Total_verify_matrix,1)
      if Total_verify_matrix(i,9)==Total_verify_matrix(i,17)
           ID_samecs_0203verify=[ID_samecs_0203verify;i];         
      end
end
Total_verify_matrix_final=Total_verify_matrix(ID_samecs_0203verify,:);
%%%%%%%%%%%%%%%
ID_training_verify=[];
for i=1:length(Training_matrix(:,16))
       Index=find(Total_verify_matrix_final(:,16)==Training_matrix(i,16));
       ID_training_verify=[ID_training_verify;Index];
end
Testing_verify_matrix_final=Total_verify_matrix_final;
Testing_verify_matrix_final(ID_training_verify,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%
Testing_verify_matrix_final(:,7);
Th=0.99;
P_th02=Testing_verify_matrix_final(:,15)>=Th;
P_th03=Testing_verify_matrix_final(:,23)>=Th;
P_th0203=P_th02+P_th03==2;
Midmatrix=Testing_verify_matrix_final;
clear Testing_verify_matrix_final
Testing_verify_matrix_final=Midmatrix(P_th0203,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_verify_ATAR0203=0;
Num_verify_AT0203=0;
num_nodetect03=0;
Id_verifynotdetect=[];
Id_verifynotdetectAT=[];
for i=1:size(Testing_verify_matrix_final,1)
    Index02=Testing_verify_matrix_final(i,16);
    Index03=Testing_verify_matrix_final(i,24);
    
    XIC02=monoXICs02_final(:,Index02);
    intervalid02=groundtruthinterval02_final(Index02,2);
    intervalmatrix02=IntervalList02_final{Index02}.intervallist;
    Ground_interval02=intervalmatrix02(intervalid02,:);
    ElutionPeak02=XIC02(Ground_interval02(1):Ground_interval02(2));

    XIC03=monoXICs03_final(:,Index03);
    intervalid03=groundtruthinterval03_final(Index03,2);
    intervalmatrix03=IntervalList03_final{Index03}.intervallist;

    for j=1:size(intervalmatrix03,1)                    
        Scanstart03=intervalmatrix03(j,1);
        Scanend03=intervalmatrix03(j,2);
        if Scanend03~=0 && Scanstart03~=0
                ElutionPeak03=XIC03(Scanstart03:Scanend03);
                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling02, newdata_sampling03, judge]=resampleforhalf(ElutionPeak02', ElutionPeak03');
                [B,BINT,R,RINT,STATS]=regress(newdata_sampling02', [ones(length(newdata_sampling03),1),newdata_sampling03']);
                Testing_verify_Rstats0203(i).AR0203(j)=STATS(1);

                Testing_verify_Rstats0203(i).AR0203_corrProb(j)=gampdf(1-STATS(1),AR0203_PARMHAT1(1),AR0203_PARMHAT1(2));
                Testing_verify_Rstats0203(i).AR0203_noncorrProb(j)=gampdf(STATS(1),AR0203_PARMHAT2(1),AR0203_PARMHAT2(2));                   

                t1=(retentiont02l1(Ground_interval02(1))+retentiont02l1(Ground_interval02(2)))/2;
                t2=(retentiont03l1(Scanstart03)+retentiont03l1(Scanend03))/2;
                t=polyval(PP23,t1)-t2;
                Testing_verify_Time0203(i).AT0203(j)=t;
                Testing_verify_Time0203(i).AT0203_corrProb(j)=normpdf(t,timemodel0203_MU_sig,timemodel0203_SIGMA_sig);
                Testing_verify_Time0203(i).AT0203_noncorrProb(j)=normpdf(t,timemodel0203_MU_notsig,timemodel0203_SIGMA_notsig);                   
        else
                num_nodetect03=num_nodetect03+1;
                Testing_verify_Rstats0203(i).AR0203(j)=0;
                Testing_verify_Rstats0203(i).AR0203_corrProb(j)=0;
                Testing_verify_Rstats0203(i).AR0203_noncorrProb(j)=0.1;                   

                Testing_verify_Time0203(i).AT0203(j)=1000000;
                Testing_verify_Time0203(i).AT0203_corrProb(j)=0;
                Testing_verify_Time0203(i).AT0203_noncorrProb(j)=0.1;                  

        end
    end
    
    [V,P]=max(log(Testing_verify_Rstats0203(i).AR0203_corrProb)+log(Testing_verify_Time0203(i).AT0203_corrProb));
%     [V,P]=max(log(Testing_verify_Rstats0203(i).AR0203_corrProb./Testing_verify_Rstats0203(i).AR0203_noncorrProb)+log(Testing_verify_Time0203(i).AT0203_corrProb./Testing_verify_Time0203(i).AT0203_noncorrProb));
    [VT,PT]=max(log(Testing_verify_Time0203(i).AT0203_corrProb));

    if P==intervalid03
        Num_verify_ATAR0203=Num_verify_ATAR0203+1;
    else Id_verifynotdetect=[Id_verifynotdetect;i];
    end
    if PT==intervalid03
        Num_verify_AT0203=Num_verify_AT0203+1;
    else Id_verifynotdetectAT=[Id_verifynotdetectAT;i];
    end
    
end
Num_verify_ATAR0203/size(Testing_verify_matrix_final,1)
Num_verify_AT0203/size(Testing_verify_matrix_final,1)

clear     Testing_verify_Rstats0203   Testing_verify_Time0203
% for oo=1:length(Id_verifynotdetect);
%     i=Id_verifynotdetect(oo);
%     Index03=Testing_verify_matrix_final(i,16);
%     intervalid03=groundtruthinterval03_final(Index03,2);
%     Testing_verify_Rstats0203(i).AR0203_corrProb
%     Testing_verify_Time0203(i).AT0203_corrProb 
%     [V,P]=max(log(Testing_verify_Rstats0203(i).AR0203_corrProb)+log(Testing_verify_Time0203(i).AT0203_corrProb));
% end




%%%%%%%%%%%%%%%%%%%%%%%%%% verify the SCIFA algorithm
ID_TF01=Totalinformationmatrix(:,8)~=0;
ID_TF02=Totalinformationmatrix(:,16)~=0;
ID_TF0102_verify=find(ID_TF01+ID_TF02==2);
length(ID_TF0102_verify)

Total_verify_matrix=Totalinformationmatrix(ID_TF0102_verify,:);
%%%%%%%%%%%% cs01 mass01 timepoint01 remarks01(from data X) intervalstart01
%%%%%%%%%%%% intervalend01 Prob01 indexindata01 (1:8)
%%%%%%%%%%%% cs01 mass01 timepoint01 remarks01(from data X) intervalstart01
%%%%%%%%%%%% intervalend01 Prob01 indexindata01 (9:16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% peplist01_final      3891                              peplist01_final  5280
%%% chargestate01_final                                  chargestate01_final
%%% groundtruthinterval01_final                     groundtruthinterval01_final
%%% IntervalList01_final                                   IntervalList01_final
%%% ms2time01_final                                         ms2time01_final
%%% monoXICs01_final                                     monoXICs01_final
%%% iso1stXICs01_final                                     iso1stXICs01_final
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% find same charge state pep in 01 01 02
ID_samecs_0102verify=[];
for i=1:size(Total_verify_matrix,1)
      if Total_verify_matrix(i,1)==Total_verify_matrix(i,9)
           ID_samecs_0102verify=[ID_samecs_0102verify;i];         
      end
end
Total_verify_matrix_final=Total_verify_matrix(ID_samecs_0102verify,:);
% %%%%%%%%%%%%%%%
% ID_training_verify=[];
% for i=1:length(Training_matrix(:,8))
%        Index=find(Total_verify_matrix_final(:,8)==Training_matrix(i,8));
%        ID_training_verify=[ID_training_verify;Index];
% end
% Testing_verify_matrix_final=Total_verify_matrix_final;
% Testing_verify_matrix_final(ID_training_verify,:)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% generate model
training_id=[];
training_timepair=[];
for i=1:size(Total_verify_matrix_final,1)
    Index01=Total_verify_matrix_final(i,8);
    Index02=Total_verify_matrix_final(i,16);
    
    XIC01=monoXICs01_final(:,Index01);
    intervalid01=groundtruthinterval01_final(Index01,2);
    intervalmatrix01=IntervalList01_final{Index01}.intervallist;
    Ground_interval01=intervalmatrix01(intervalid01,:);
    ElutionPeak01=XIC01(Ground_interval01(1):Ground_interval01(2));

    XIC02=monoXICs02_final(:,Index02);
    intervalid02=groundtruthinterval02_final(Index02,2);
    intervalmatrix02=IntervalList02_final{Index02}.intervallist;
    Ground_interval02=intervalmatrix02(intervalid02,:);
    ElutionPeak02=XIC02(Ground_interval02(1):Ground_interval02(2));

    
    Int01_max=max(XIC01(Ground_interval01(1):Ground_interval01(2)));
    Int02_max=max(XIC02(Ground_interval02(1):Ground_interval02(2)));

    if Int01_max>=10000000 && Int02_max>=10000000
       training_id=[training_id; i];
       training_timepair=[training_timepair;Total_verify_matrix_final(i,3),Total_verify_matrix_final(i,11)];
    end
end
length(training_id)
Training_verify_matrix_final=Total_verify_matrix_final(training_id,:);

PP=polyfit(training_timepair(:,1),training_timepair(:,2),4);
ARdiff_corr=[];ARdiff_non_corr=[];
ATdiff_corr=[];ATdiff_non_corr=[];

for i=1:size(Training_verify_matrix_final,1)
    
    Index01=Total_verify_matrix_final(i,8);
    Index02=Total_verify_matrix_final(i,16);
    
    XIC01=monoXICs01_final(:,Index01);
    intervalid01=groundtruthinterval01_final(Index01,2);
    intervalmatrix01=IntervalList01_final{Index01}.intervallist;
    Ground_interval01=intervalmatrix01(intervalid01,:);
    ElutionPeak01=XIC01(Ground_interval01(1):Ground_interval01(2));

    XIC02=monoXICs02_final(:,Index02);
    intervalid02=groundtruthinterval02_final(Index02,2);
    intervalmatrix02=IntervalList02_final{Index02}.intervallist;
    Ground_interval02=intervalmatrix02(intervalid02,:);
    ElutionPeak02=XIC02(Ground_interval02(1):Ground_interval02(2));  

    %%%%%% AR AT 0102
    for j1=1:size(intervalmatrix01,1)
        for j2=1:size(intervalmatrix02,1)

                scanstart01=intervalmatrix01(j1,1);
                scanend01=intervalmatrix01(j1,2);
                monoelutionprofile01=XIC01(scanstart01:scanend01);

                scanstart02=intervalmatrix02(j2,1);
                scanend02=intervalmatrix02(j2,2);
                monoelutionprofile02=XIC02(scanstart02:scanend02);

                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf( monoelutionprofile01', monoelutionprofile02');
                [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                Training_Rstats0102(i).AR0102(j1,j2)=STATS(1);

                t1=(retentiont01l1(scanstart01)+retentiont01l1(scanend01))/2;
                t2=(retentiont02l1(scanstart02)+retentiont02l1(scanend02))/2;
                t=polyval(PP,t1)-t2;
                Training_Time0102(i).AT0102(j1,j2)=t;

        end
    end
    %%%%%% 
    
    ARdiff_corr=[ARdiff_corr,Training_Rstats0102(i).AR0102(intervalid01,intervalid02)];
    ATdiff_corr=[ATdiff_corr,Training_Time0102(i).AT0102(intervalid01,intervalid02)];

    k1=intervalid01;
    k2=intervalid02;
    tdiff=Training_Time0102(i).AT0102(k1,:);
    ardiff=Training_Rstats0102(i).AR0102(k1,:);
    tdiff(k2)=[];      ardiff(k2)=[];
    tdiffv1=reshape(tdiff,1,length(tdiff));
    ATdiff_non_corr=[ATdiff_non_corr,tdiffv1];
    ARdiff_non_corr=[ARdiff_non_corr,ardiff]; 
    
end
%%%%%% guassian fit timeshift
[timemodel_MU_sig,timemodel_SIGMA_sig] = normfit(ATdiff_corr);
[timemodel_MU_notsig,timemodel_SIGMA_notsig] = normfit(ATdiff_non_corr);
%%%%%% gamma fit AR
aa=find(ARdiff_non_corr>=0 & ARdiff_non_corr<=1);
ARdiff_non_corr_noNaN=ARdiff_non_corr(aa);
aa=find(ARdiff_corr>=0 & ARdiff_corr<=1);
ARdiff_corr_noNaN=ARdiff_corr(aa);
AR_PARMHAT1=gamfit(1-ARdiff_corr_noNaN);
AR_PARMHAT2=gamfit(ARdiff_non_corr_noNaN);

%%%% model_parameters:
%%%% PP timemodel_MU_sig timemodel_SIGMA_sig
%%%% timemodel_MU_notsig timemodel_SIGMA_notsig
%%%% AR_PARMHAT1 AR_PARMHAT2

%%%%%%%%%
Testing_verify_matrix_final=Total_verify_matrix_final;
Testing_verify_matrix_final(training_id,:)=[];
Testing_verify_matrix_final(:,7);
Th=0.99;
P_th01=Testing_verify_matrix_final(:,7)>=Th;
P_th02=Testing_verify_matrix_final(:,15)>=Th;
P_th0102=P_th01+P_th02==2;
Midmatrix=Testing_verify_matrix_final;
clear Testing_verify_matrix_final
Testing_verify_matrix_final=Midmatrix(P_th0102,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_verify_ATAR0102=0;
Num_verify_AT0102=0;
num_nodetect02=0;
Id_verifynotdetect=[];
Id_verifynotdetectAT=[];
for i=1:size(Testing_verify_matrix_final,1)
    Index01=Testing_verify_matrix_final(i,8);
    Index02=Testing_verify_matrix_final(i,16);
    
    XIC01=monoXICs01_final(:,Index01);
    intervalid01=groundtruthinterval01_final(Index01,2);
    intervalmatrix01=IntervalList01_final{Index01}.intervallist;
    Ground_interval01=intervalmatrix01(intervalid01,:);
    ElutionPeak01=XIC01(Ground_interval01(1):Ground_interval01(2));

    XIC02=monoXICs02_final(:,Index02);
    intervalid02=groundtruthinterval02_final(Index02,2);
    intervalmatrix02=IntervalList02_final{Index02}.intervallist;

    for j=1:size(intervalmatrix02,1)                    
        Scanstart02=intervalmatrix02(j,1);
        Scanend02=intervalmatrix02(j,2);
        if Scanend02~=0 && Scanstart02~=0
                ElutionPeak02=XIC02(Scanstart02:Scanend02);
                [newdata_nosamp, newdata_sampwithoutshift, newdata_sampling01, newdata_sampling02, judge]=resampleforhalf(ElutionPeak01', ElutionPeak02');
                [B,BINT,R,RINT,STATS]=regress(newdata_sampling01', [ones(length(newdata_sampling02),1),newdata_sampling02']);
                Testing_verify_Rstats0102(i).AR0102(j)=STATS(1);

                Testing_verify_Rstats0102(i).AR0102_corrProb(j)=gampdf(1-STATS(1),AR_PARMHAT1(1),AR_PARMHAT1(2));
                Testing_verify_Rstats0102(i).AR0102_noncorrProb(j)=gampdf(STATS(1),AR_PARMHAT2(1),AR_PARMHAT2(2));                   

                t1=(retentiont01l1(Ground_interval01(1))+retentiont01l1(Ground_interval01(2)))/2;
                t2=(retentiont02l1(Scanstart02)+retentiont02l1(Scanend02))/2;
                t=polyval(PP,t1)-t2;
                Testing_verify_Time0102(i).AT0102(j)=t;
                Testing_verify_Time0102(i).AT0102_corrProb(j)=normpdf(t,timemodel_MU_sig,timemodel_SIGMA_sig);
                Testing_verify_Time0102(i).AT0102_noncorrProb(j)=normpdf(t,timemodel_MU_notsig,timemodel_SIGMA_notsig);                   
        else
                num_nodetect02=num_nodetect02+1;
                Testing_verify_Rstats0102(i).AR0102(j)=0;
                Testing_verify_Rstats0102(i).AR0102_corrProb(j)=0;
                Testing_verify_Rstats0102(i).AR0102_noncorrProb(j)=0.1;                   

                Testing_verify_Time0102(i).AT0102(j)=1000000;
                Testing_verify_Time0102(i).AT0102_corrProb(j)=0;
                Testing_verify_Time0102(i).AT0102_noncorrProb(j)=0.1;                  

        end
    end
    
    [V,P]=max(log(Testing_verify_Rstats0102(i).AR0102_corrProb)+log(Testing_verify_Time0102(i).AT0102_corrProb));
%     [V,P]=max(log(Testing_verify_Rstats0102(i).AR0102_corrProb./Testing_verify_Rstats0102(i).AR0102_noncorrProb)+log(Testing_verify_Time0102(i).AT0102_corrProb./Testing_verify_Time0102(i).AT0102_noncorrProb));
    [VT,PT]=max(log(Testing_verify_Time0102(i).AT0102_corrProb));

    if P==intervalid02
        Num_verify_ATAR0102=Num_verify_ATAR0102+1;
    else Id_verifynotdetect=[Id_verifynotdetect;i];
    end
    if PT==intervalid02
        Num_verify_AT0102=Num_verify_AT0102+1;
    else Id_verifynotdetectAT=[Id_verifynotdetectAT;i];
    end
    
end

Num_verify_ATAR0102/size(Testing_verify_matrix_final,1)
Num_verify_AT0102/size(Testing_verify_matrix_final,1)

clear     Testing_verify_Rstats0102   Testing_verify_Time0102




