clc
clear all



save_path='C:\Users\Long\Desktop\SILAC Processing\Result2';
%%%%%%%%%%%%%%%%%% read Maxquant txt file and generate Information matrix
%%%%%%%%%%%%%%% Total peptide information
%%%%%%%%%%%%%%% column name: peptide_seq cs mass  ms2scannumber HLR protein
%%%%%%%%%%%%%%% totalmzlist iso (1-8 column for 1st replicate 9-16 for 2nd 17-24 for 3rd)
%%%%%% note: load Maxquant processing result of Pre-ratio replicates (3 data with same HLR)
%%%%%% the sum of iso is not normalized to 1.

path='C:\Users\Long\Desktop\SILAC Processing\CC35_data\combined\txt';
[Information_Matrix_K_labeled_A, Information_Matrix_A]=readMaxquantResult(path);
path='C:\Users\Long\Desktop\SILAC Processing\CC50_data\combined\txt';
[Information_Matrix_K_labeled_B, Information_Matrix_B]=readMaxquantResult(path);
path='C:\Users\Long\Desktop\SILAC Processing\CC50_data\combined\txt';
[Information_Matrix_K_labeled_C, Information_Matrix_C]=readMaxquantResult(path);

%%%%%%%%%% combine information matrix from 3 input matrix (Remove redundant entries)
Information_Matrix_Total=CombineOriginalInformationMatrix(Information_Matrix_K_labeled_A,Information_Matrix_K_labeled_B);
Information_Matrix_Total=CombineOriginalInformationMatrix(Information_Matrix_Total,Information_Matrix_K_labeled_C);
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%% read mzxml data 
% MzxmlFile_path='C:\Users\Long\Desktop\SILAC Processing\mzXML\CC35.mzXML';
% [retentiontl1_A,datal1_A,peakl_A,retentiont_A,MZInt_l1l2_A]=readrawdata(MzxmlFile_path);
% save([save_path,'\Mzxml_A_readin'], 'retentiontl1_A','peakl_A','retentiont_A');
% MzxmlFile_path='C:\Users\Long\Desktop\SILAC Processing\mzXML\CC50.mzXML';
% [retentiontl1_B,datal1_B,peakl_B,retentiont_B,MZInt_l1l2_B]=readrawdata(MzxmlFile_path);
% save([save_path,'\Mzxml_B_readin'], 'retentiontl1_B','peakl_B','retentiont_B');
% MzxmlFile_path='C:\Users\Long\Desktop\SILAC Processing\mzXML\CC70.mzXML';
% [retentiontl1_C,datal1_C,peakl_C,retentiont_C,MZInt_l1l2_C]=readrawdata(MzxmlFile_path);
% save([save_path,'\Mzxml_C_readin'], 'retentiontl1_C','peakl_C','retentiont_C');

load([save_path,'\Mzxml_A_readin']);
load([save_path,'\Mzxml_B_readin']);
load([save_path,'\Mzxml_C_readin']);

TOF_totalmzList=zeros(size(Information_Matrix_Total,1)*8,1);
TOTAL_isolist=zeros(size(Information_Matrix_Total,1),4);
TOTAL_cs=zeros(size(Information_Matrix_Total,1),1);
for i=1:size(Information_Matrix_Total,1)
    Massseq=[Information_Matrix_Total{i,3},Information_Matrix_Total{i,11},Information_Matrix_Total{i,19}];
    ID_n0=find(Massseq~=0);
    TOF_totalmzList((i-1)*8+1:i*8)=Information_Matrix_Total{i,(ID_n0(1)-1)*8+7};
    TOTAL_isolist(i,:)=Information_Matrix_Total{i,(ID_n0(1)-1)*8+8};
    TOTAL_cs(i)=Information_Matrix_Total{i,(ID_n0(1)-1)*8+2}(1);
end
%%%%%%%generate XICs (The ith peptide has 8 XICs continuous in colomun IDs (i-1)*8+1:i*8)
tolerance=20;
XICs_A=getXIC_LC_new(peakl_A,TOF_totalmzList,tolerance);
XICs_B=getXIC_LC_new(peakl_B,TOF_totalmzList,tolerance);
XICs_C=getXIC_LC_new(peakl_C,TOF_totalmzList,tolerance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% interval detection
%%%%%%%%%%%%%%%%% the returned data structure contains intervals of LC peaks. the
%%%%%%%%%%%%%%%%% final intervallist is called intervallist_after7_combine
pep_number=size(XICs_A,2)/8;
for i=1:pep_number
    isoList=TOTAL_isolist(i,:);
    cs=TOTAL_cs(i);
    for filecount=1:3
       if filecount==1 
           TOF_XICs=XICs_A;
           Retent=retentiontl1_A;
       elseif  filecount==2
            TOF_XICs=XICs_B;
            Retent=retentiontl1_B;
       elseif filecount==3
            TOF_XICs=XICs_C;
            Retent=retentiontl1_C;
       end
        pepXICs=TOF_XICs(:,(i-1)*8+1:i*8);
        [IntervalList,Jud_detect_good]=SILAC_Verification_intervaldetection_Spec_Pepv1(isoList,pepXICs);
        IntervalListBig{i,filecount}=IntervalList;
        IntervalListBig{i,filecount}.chargeState=cs;
         
%         if IntervalList.totalInterval>0    
%             colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
%             figure
%             
%             for index=1:6
%                 plot(Retent,pepXICs(:,index),colorarray(index))
%                 hold on;
%             end
%             height=1000;
%             stem(T_ms2,height*2)
%             stem(Retent(IntervalList.intervallist_after7_combine(:,1)),2*height*ones(length(IntervalList.intervallist_after7_combine(:,1)),1),'ro')
%             stem(Retent(IntervalList.intervallist_after7_combine(:,2)),2*height*ones(length(IntervalList.intervallist_after7_combine(:,2)),1),'ko')
% %             stem(retentiontl1_A(IntervalList.intervallist_after7(:,1)),1.5*height*ones(length(IntervalList.intervallist_after7(:,1)),1),'r*')
% %             stem(retentiontl1_A(IntervalList.intervallist_after7(:,2)),1.5*height*ones(length(IntervalList.intervallist_after7(:,2)),1),'k*')
% %             stem(retentiontl1_A(IntervalList.intervallist_after6(:,1)),1*height*ones(length(IntervalList.intervallist_after6(:,1)),1),'rs')
% %             stem(retentiontl1_A(IntervalList.intervallist_after6(:,2)),1*height*ones(length(IntervalList.intervallist_after6(:,2)),1),'ms')
%         end
%         intervalNumbers(i,filecount)=IntervalList.totalInterval;
    end     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find training data set
%%%%%%%%% 1. find intersection of the peptide with ms2 time point
%%%%%%%%%% Information_Matrix_Total column name: 1[peptide_seq] 2[cs]
%%%%%%%%%% 3[mass] 4[ms2scannumber] 5[HLR] 6[protein] 7[totalmzlist] 8[iso]
%%%%%%%%%% INTERVAL_Matrix column name: 1[filename] 2[interval id] 3[scanstart] 4[scanend]
%%%%%%%%%% each peptide has an entry in INTERVAL_Matrix.
%%%%%%%%%% file name indicates in which replicate the peptide is picked up
%%%%%%%%%% by ms2 the value is 1/0,2/0 or 3/0; interval id is the index of
%%%%%%%%%% the ms2 LC peak in the interval list.
%%%%%%%%%% the out put is the INTERVAL_Matrix.
INTERVAL_Matrix=zeros(size(Information_Matrix_Total,1),12);
DetectNumberinA=0; DetectNumberinB=0; DetectNumberinC=0;
for i=1:pep_number
        MS2Scannumber_Vector=[Information_Matrix_Total{i,4},Information_Matrix_Total{i,12},Information_Matrix_Total{i,20}];
        ID_n0=find(MS2Scannumber_Vector~=0);
        for j=1:length(ID_n0) 
                filecount=ID_n0(j);
                if filecount==1
                    retentiontl1=retentiontl1_A;
                    retentiont=retentiont_A;
                else if filecount==2
                            retentiontl1=retentiontl1_B;
                            retentiont=retentiont_B;
                    else if filecount==3                            
                                retentiontl1=retentiontl1_C;
                                retentiont=retentiont_C;
                        end
                    end
                end  
                
                MS2Scan=MS2Scannumber_Vector(filecount);
                Intervals=IntervalListBig{i,filecount}.intervallist_after7_combine;  
                if sum(sum(Intervals))~=0
                    retentiontl1v1=retentiontl1';
                    Pep_ms2time=retentiont(MS2Scan,1);
                    INTD=retentiontl1v1(Intervals)-Pep_ms2time;
                    ID_MS2=find(INTD(:,1).*INTD(:,2)<=0);
                    if ~isempty(ID_MS2)
                        INTERVAL_Matrix(i,(filecount-1)*4+1)=filecount;
                        INTERVAL_Matrix(i,(filecount-1)*4+2)=ID_MS2;
                        INTERVAL_Matrix(i,(filecount-1)*4+3)=Intervals(ID_MS2,1);
                        INTERVAL_Matrix(i,(filecount-1)*4+4)=Intervals(ID_MS2,2);
                        if filecount==1
                            DetectNumberinA=DetectNumberinA+1; 
                        else if filecount==2
                            DetectNumberinB=DetectNumberinB+1; 
                            else if filecount==3
                            DetectNumberinC=DetectNumberinC+1; 
                                end
                            end
                        end
                    end            
                end
        end
end

% %%%%%%%%%%% calculate the intersection of the A B C ms2 peptides datasets
% Count=sum([INTERVAL_Matrix(:,1),INTERVAL_Matrix(:,5)/2,INTERVAL_Matrix(:,9)/3],2);
% CountAB=sum([INTERVAL_Matrix(:,1),INTERVAL_Matrix(:,5)/2],2);
% CountBC=sum([INTERVAL_Matrix(:,5)/2,INTERVAL_Matrix(:,9)/3],2);
% CountAC=sum([INTERVAL_Matrix(:,1),INTERVAL_Matrix(:,9)/3],2);
% IDwithms2interval=find(Count~=0);
% IDwith3ms2interval=find(Count==3);
% IDwithms2intervalAB=find(CountAB==2);
% IDwithms2intervalBC=find(CountBC==2);
% IDwithms2intervalAC=find(CountAC==2);
% DetectNumberinA
% DetectNumberinB
% DetectNumberinC
% DetectNumberUnion=length(IDwithms2interval)
% DetectNumberintersection=length(IDwith3ms2interval)
% DetectNumberintersectionAB=length(IDwithms2intervalAB)
% DetectNumberintersectionBC=length(IDwithms2intervalBC)
% DetectNumberintersectionAC=length(IDwithms2intervalAC)

%%%%%%%%%%%%%%%%%%%Check interval
%%%%%%%%%%%%%%%%%%% diagnostic code for interval detection (can be commented out) 
% for i=1:pep_number
%     IDVectorwithMS2=[INTERVAL_Matrix(i,1),INTERVAL_Matrix(i,5),INTERVAL_Matrix(i,9)];
%     ms2ScanVector=[Information_Matrix_Total{i,4},Information_Matrix_Total{i,12},Information_Matrix_Total{i,20}];
%     for filecount=1:3
%        if filecount==1 
%            TOF_XICs=XICs_A;
%            Retent=retentiontl1_A;
%            Retentms2=retentiont_A;
%        elseif  filecount==2
%             TOF_XICs=XICs_B;
%             Retent=retentiontl1_B;
%             Retentms2=retentiont_B;
%        elseif filecount==3
%             TOF_XICs=XICs_C;
%             Retent=retentiontl1_C;
%             Retentms2=retentiont_C;
%        end
%         pepXICs=TOF_XICs(:,(i-1)*8+1:i*8);
%         IntervalList=IntervalListBig{i,filecount};
%        
%         colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
%         figure            
%         for index=1:6
%             plot(Retent,pepXICs(:,index),colorarray(index))
%             hold on;
%         end
%         height=10000;
%         if ms2ScanVector(filecount)~=0
%             T_ms2=Retentms2(ms2ScanVector(filecount));
%             Interval_M=IntervalList.intervallist_after7_combine;
%         else
%             T_ms2=Retentms2(ms2ScanVector(filecount)+1);
%             Interval_M=IntervalList.intervallist_after7_combine+1;
%         end
%         stem(T_ms2,height*2,'r*')
%         stem(Retent(Interval_M(:,1)),2*height*ones(length(Interval_M(:,1)),1),'ro')
%         stem(Retent(Interval_M(:,2)),2*height*ones(length(Interval_M(:,2)),1),'ko')
% %             stem(retentiontl1_A(IntervalList.intervallist_after7(:,1)),1.5*height*ones(length(IntervalList.intervallist_after7(:,1)),1),'r*')
% %             stem(retentiontl1_A(IntervalList.intervallist_after7(:,2)),1.5*height*ones(length(IntervalList.intervallist_after7(:,2)),1),'k*')
% %             stem(retentiontl1_A(IntervalList.intervallist_after6(:,1)),1*height*ones(length(IntervalList.intervallist_after6(:,1)),1),'rs')
% %             stem(retentiontl1_A(IntervalList.intervallist_after6(:,2)),1*height*ones(length(IntervalList.intervallist_after6(:,2)),1),'ms')
%     end
%     
% end
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%need to add to detection wheather the XICs and interval
%%%%%%%%%%%%%information exist or not
save([save_path,'\XICandIntervalInformation'], 'XICs_A', 'XICs_B', 'XICs_C', 'IntervalListBig', 'INTERVAL_Matrix', 'Information_Matrix_Total')
load([save_path,'\XICandIntervalInformation']);
%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% select training data (peptides with intervals that contains 3 ms2 time points)
IDMwithMS2=[INTERVAL_Matrix(:,1)/1,INTERVAL_Matrix(:,5)/2,INTERVAL_Matrix(:,9)/3];
IDVwithMS2=sum(IDMwithMS2,2);
Training_ID=find(IDVwithMS2==3);
%%%%%%%%%%% build models AT AR AKL (AB AC BC datasets) and save them in
%%%%%%%%%%% Model_Training with AB AC BC order

%%%%%%%%%%%%%%%%%%%%%%%Generate Warping function of AB AC BC
for i=1:length(Training_ID)
    P_ID=Training_ID(i);
    MS2scan_number_vector=[Information_Matrix_Total{P_ID,4},Information_Matrix_Total{P_ID,12},Information_Matrix_Total{P_ID,20}];
    %%%%%%%%%%%%%%%% combine peptide information of A B C before building
    %%%%%%%%%%%%%%%% models
    for filecount=1:3
        if filecount==1 
           Retent=retentiontl1_A;
           Retentms2=retentiont_A;
       elseif  filecount==2
            Retent=retentiontl1_B;
            Retentms2=retentiont_B;
       elseif filecount==3
            Retent=retentiontl1_C;
            Retentms2=retentiont_C;
        end
        Interval_after7_combine=IntervalListBig{P_ID,filecount}.intervallist_after7_combine;
        MS2_ID=INTERVAL_Matrix(P_ID,(filecount-1)*4+2);
        MS2_SC=MS2scan_number_vector(filecount);
        
        MS2_Time_Matrix(i,filecount)=Retentms2(MS2_SC);
        MS2_Time_Matrix_Scan(i,filecount)=sum(Retent(Interval_after7_combine(MS2_ID,:)))/2;
    end   
end

PP_Parameter.PP_AB=polyfit(MS2_Time_Matrix(:,1),MS2_Time_Matrix(:,2),4);
PP_Parameter.PP_AC=polyfit(MS2_Time_Matrix(:,1),MS2_Time_Matrix(:,3),4);
PP_Parameter.PP_BC=polyfit(MS2_Time_Matrix(:,2),MS2_Time_Matrix(:,3),4);

% X=3800:10:5400;
% Y=polyval(PP_AB,X);
% figure
% plot(MS2_Time_Matrix(:,1),MS2_Time_Matrix(:,2),'r*');hold on
% plot(X,Y,'k');hold on
% plot(MS2_Time_Matrix_Scan(:,1),MS2_Time_Matrix_Scan(:,2),'ko')

% find(abs(MS2_Time_Matrix(:,2)-MS2_Time_Matrix_Scan(:,2))>=200)
%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Training_ID)
    P_ID=Training_ID(i);
    IsoList=Information_Matrix_Total{P_ID,8};
    MS2scan_number_vector=[Information_Matrix_Total{P_ID,4},Information_Matrix_Total{P_ID,12},Information_Matrix_Total{P_ID,20}];
    %%%%%%%%%%%%%%%% combine peptide information of A B C before building
    %%%%%%%%%%%%%%%% models
    for filecount=1:3
        if filecount==1 
           TOF_XICs=XICs_A;
           Retent=retentiontl1_A;
           Retentms2=retentiont_A;
       elseif  filecount==2
            TOF_XICs=XICs_B;
            Retent=retentiontl1_B;
            Retentms2=retentiont_B;
       elseif filecount==3
            TOF_XICs=XICs_C;
            Retent=retentiontl1_C;
            Retentms2=retentiont_C;
        end
        
       Interval_after7_combine=IntervalListBig{P_ID,filecount}.intervallist_after7_combine;
       MS2_ID=INTERVAL_Matrix(P_ID,(filecount-1)*4+2);
       MS2_SC=MS2scan_number_vector(filecount);        
        
       for Interval_Row=1:size(Interval_after7_combine,1)
           IntervalElutionMatrix{Interval_Row}=TOF_XICs(Interval_after7_combine(Interval_Row,1):Interval_after7_combine(Interval_Row,2),(P_ID-1)*8+1:P_ID*8);
           Interval_after7_combine_Time{Interval_Row}=Retent(Interval_after7_combine(Interval_Row,1):Interval_after7_combine(Interval_Row,2));
       end
        
        Training_Information_Matrix{i,filecount}.iso=IsoList;
        Training_Information_Matrix{i,filecount}.Interval_after7_combine=Interval_after7_combine;
        Training_Information_Matrix{i,filecount}.MS2_ID=MS2_ID;
        Training_Information_Matrix{i,filecount}.Interval_after7_combine_Time=Interval_after7_combine_Time;
        Training_Information_Matrix{i,filecount}.MS2timePoint=Retentms2(MS2_SC);
        Training_Information_Matrix{i,filecount}.IntervalElutionMatrix=IntervalElutionMatrix;        
    end

    TrainingInformation=Training_Information_Matrix(i,:);
    Model_Training{i}=Generate_Training_Scores(TrainingInformation,PP_Parameter);
end

save([save_path,'\Training_Information_Matrix'], 'Training_Information_Matrix')
save([save_path,'\Model_Training'], 'Model_Training')


Parameters=BuildATARAKLModelParameters(Training_Information_Matrix,Model_Training);

save([save_path,'\Model_Training'], 'PP_Parameter')
save([save_path,'\Parameters'], 'Parameters')
load([save_path,'\Model_Training'], 'PP_Parameter')
load([save_path,'\Parameters'], 'Parameters')
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%Generate Testing datasets

%%%%% build time model and score seperately from the AR AKL 
%%%%% The warping parameter are A->B A->C and B->C
%%%%% Based on the MS2scan_number_vector, we should get two pairs F1<->B1
%%%%% and F2<->B2.
for i=1:size(Information_Matrix_Total,1)
    IsoList=Information_Matrix_Total{i,8};
    ID_notzeros=[Information_Matrix_Total{i,4},Information_Matrix_Total{i,12},Information_Matrix_Total{i,20}]~=0;
    MS2scan_number_vector=zeros(1,3);
    MS2scan_number_vector(ID_notzeros)=1;    
    
     for filecount=1:3
        if filecount==1 
           TOF_XICs=XICs_A;
           Retent=retentiontl1_A;
           Retentms2=retentiont_A;
       elseif  filecount==2
            TOF_XICs=XICs_B;
            Retent=retentiontl1_B;
            Retentms2=retentiont_B;
       elseif filecount==3
            TOF_XICs=XICs_C;
            Retent=retentiontl1_C;
            Retentms2=retentiont_C;
       end
        
       Interval_after7_combine=IntervalListBig{i,filecount}.intervallist_after7_combine;
       MS2_ID=INTERVAL_Matrix(i,(filecount-1)*4+2);
       MS2_SC=MS2scan_number_vector(filecount);
       if sum(sum(Interval_after7_combine))==0
           IntervalElutionMatrix{1}=0;
           Interval_after7_combine_Time=0;
           Scan_max=0;
           Time_max_point=0;
       else
           for Interval_Row=1:size(Interval_after7_combine,1)
               IntervalElutionMatrix{Interval_Row}=TOF_XICs(Interval_after7_combine(Interval_Row,1):Interval_after7_combine(Interval_Row,2),(i-1)*8+1:i*8);
               Interval_after7_combine_Time{Interval_Row}=Retent(Interval_after7_combine(Interval_Row,1):Interval_after7_combine(Interval_Row,2));
               [a,b]=max(IntervalElutionMatrix{Interval_Row});
               [c,d]=max(a);
               Scan_max(Interval_Row)=b(d);
               Time_max_point(Interval_Row)=Interval_after7_combine_Time{Interval_Row}(b(d));
           end       
       end
        
        Total_test_Information_Matrix{i,filecount}.iso=IsoList;
        Total_test_Information_Matrix{i,filecount}.Interval_after7_combine=Interval_after7_combine;
        Total_test_Information_Matrix{i,filecount}.MS2_ID=MS2_ID;
        Total_test_Information_Matrix{i,filecount}.Interval_after7_combine_Time=Interval_after7_combine_Time;
        Total_test_Information_Matrix{i,filecount}.Scan_max=Scan_max;
        Total_test_Information_Matrix{i,filecount}.Time_max_point=Time_max_point;
        if MS2_SC==0
             Total_test_Information_Matrix{i,filecount}.MS2timePoint=0;
        else
            Total_test_Information_Matrix{i,filecount}.MS2timePoint=Retentms2(MS2_SC);
        end
        Total_test_Information_Matrix{i,filecount}.IntervalElutionMatrix=IntervalElutionMatrix;
        
        clear Interval_after7_combine_Time IntervalElutionMatrix Interval_after7_combine Scan_max Time_max_point
     end
    
end

for i=1:size(Information_Matrix_Total,1)
    J_Val=[INTERVAL_Matrix(i,1),INTERVAL_Matrix(i,5),INTERVAL_Matrix(i,9)];
    J_Vec=J_Val~=0;
    J_V=J_Vec*[1;2;4];
    switch J_V
        case 0 %%[A B C]=[0 0 0]

            Row01=0;            Col01=0;            Row02=0;            Col02=0;
            RowElutionProf01{1}=0;            ColElutionProf01{1}=0;
            RowElutionProf02{1}=0;            ColElutionProf02{1}=0;
            
        case 1 %%[A B C]=[1 0 0] (Time warping PP(A)-B and PP(A)-C)
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 2 %%[A B C]=[0 1 0] (Time warping B-PP(A) and PP(B)-C)
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;  
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 3 %%[A B C]=[1 1 0] (Time warping PP(A)-C and PP(B)-C)
            Row01=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,3}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 4 %%[A B C]=[0 0 1] (Time warping C-PP(A) and C-PP(B))
            Row01=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,3}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point; 
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 5 %%[A B C]=[1 0 1] (Time warping PP(A)-B and C-PP(B))
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point; 
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 6 %%[A B C]=[0 1 1] (Time warping B-PP(A) and C-PP(A))
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 7 %%[A B C]=[1 1 1]    
            %%%%%%%%%%%%%% it is not necessary to do alignment in this
            %%%%%%%%%%%%%% situation
            Row01=0;            Col01=0;            Row02=0;            Col02=0;
            RowElutionProf01{1}=0;            ColElutionProf01{1}=0;
            RowElutionProf02{1}=0;            ColElutionProf02{1}=0;
    end
    TestModel{i}.J_Vec=J_Vec;
    TestModel{i}.Trow01=Row01;
    TestModel{i}.Trow02=Row02;
    TestModel{i}.Tcol01=Col01;
    TestModel{i}.Tcol02=Col02;
    TestModel{i}.Tmatrix01=Row01'*ones(1,length(Col01))-ones(length(Row01),1)*Col01;
    TestModel{i}.Tmatrix02=Row02'*ones(1,length(Col02))-ones(length(Row02),1)*Col02;
    
    [TestModel{i}.ARmatrix01,TestModel{i}.AKLmatrix01]=CalARAKLMatrix(RowElutionProf01,ColElutionProf01);
    [TestModel{i}.ARmatrix02,TestModel{i}.AKLmatrix02]=CalARAKLMatrix(RowElutionProf02,ColElutionProf02);
    
    for kk=1:3
        TestModel{i}.IntervalListInformation{kk}=IntervalListBig{i,kk};
        TestModel{i}.ImportantInf{kk}=Total_test_Information_Matrix{i,kk};        
    end
    
    clear RowElutionProf01 ColElutionProf01 RowElutionProf02 ColElutionProf02
end
    
for i=1:length(TestModel)
    J_Vec=TestModel{i}.J_Vec;
    J_V=J_Vec*[1;2;4];
    switch J_V
        case 0 %%[A B C]=[0 0 0]
        %%%%%%%%%%%% not the right situation that can align without ms2
        %%%%%%%%%%%% information                    
        case 1 %%[A B C]=[1 0 0] (Time warping PP(A)-B and PP(A)-C)
            Para01=Parameters{1};Para02=Parameters{2};
        case 2 %%[A B C]=[0 1 0] (Time warping B-PP(A) and PP(B)-C)
            Para01=Parameters{1};Para02=Parameters{3};
        case 3 %%[A B C]=[1 1 0] (Time warping PP(A)-C and PP(B)-C)
            Para01=Parameters{2};Para02=Parameters{3};
        case 4 %%[A B C]=[0 0 1] (Time warping C-PP(A) and C-PP(B))
            Para01=Parameters{2};Para02=Parameters{3};
        case 5 %%[A B C]=[1 0 1] (Time warping PP(A)-B and C-PP(B))
            Para01=Parameters{1};Para02=Parameters{3};
        case 6 %%[A B C]=[0 1 1] (Time warping B-PP(A) and C-PP(A))
            Para01=Parameters{1};Para02=Parameters{2};
        case 7 %%[A B C]=[1 1 1]    
            %%%%%%%%%%%%%% it is not necessary to do alignment in this
            %%%%%%%%%%%%%% situation
    end
    if J_V~=0 && J_V~=7
        
        [TestModel{i}.ATScoreMatrixCorr01,TestModel{i}.ATScoreMatrixNCorr01]=GenATScore(TestModel{i}.Tmatrix01,Para01);
        [TestModel{i}.ATScoreMatrixCorr02,TestModel{i}.ATScoreMatrixNCorr02]=GenATScore(TestModel{i}.Tmatrix02,Para02);
        [TestModel{i}.ARScoreMatrixCorr01,TestModel{i}.ARScoreMatrixNCorr01]=GenARScore(TestModel{i}.ARmatrix01,Para01);
        [TestModel{i}.ARScoreMatrixCorr02,TestModel{i}.ARScoreMatrixNCorr02]=GenARScore(TestModel{i}.ARmatrix02,Para02);
        [TestModel{i}.AKLScoreMatrixCorr01,TestModel{i}.AKLScoreMatrixNCorr01]=GenAKLScore(TestModel{i}.AKLmatrix01,Para01);
        [TestModel{i}.AKLScoreMatrixCorr02,TestModel{i}.AKLScoreMatrixNCorr02]=GenAKLScore(TestModel{i}.AKLmatrix02,Para02);        
        
    end
end
save([save_path,'\TestModel'], 'TestModel')
load([save_path,'\TestModel'])
load([save_path,'\XICandIntervalInformation'])
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% plot interval elution profile
%%%%%%%%%%%%%%%%%%% diagnostic code for interval detection (can be commented out) 
% for i=5:10%:length(TestModel)  
%     colorarray=['r', 'k', 'g', 'b', 'm', 'y', 'r+', 'k+'];
%     for filecount=1:3
%         ID_ms2=TestModel{i}.ImportantInf{filecount}.MS2_ID;
%         Elution_Matrix=TestModel{i}.ImportantInf{filecount}.IntervalElutionMatrix{ID_ms2};
% %         Elution_vector=[];
%         figure
%         for j=1:8                    
%             plot(Elution_Matrix(:,j),colorarray(j))
%             hold on
%         end
%         clear Elution_Matrix
%     end    
% end    
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% pick the aligned pairs
ID_nondetected=[];
for i=1:length(TestModel)
        FinalResult{i}.ImportantInf=TestModel{i}.ImportantInf;
        %%%%%%%%%%%%%%% ScoreMatrix 7 cell for AT AR AKL AT+AR AT+AKL
        %%%%%%%%%%%%%%% AR+AKL AT+AR+AKL
        %%%%%%%%%%%%%%% 01 for 1st data pair 02 for 2nd data pair
        
        ChooseScoreMatrix=5; %%%%% 1:AT 2:AR 3:AKL 4:AT+AR 5:AT+AKL
        %%%%%%%%%%%%%%% 6:AR+AKL 7:AT+AR+AKL
        J_Vec=TestModel{i}.J_Vec;
        J_Val=sum(J_Vec*[1 2 4]');
        switch J_Val
            
            case 0%%[0 0 0]
              ID_nondetected=[ID_nondetected;i];
              FinalResult{i}.J_Vec=J_Vec;
              peptide{i} = '-';
              protein{i} = '-';
              iso{i} = 0;
            case 1 %%[1 0 0]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});  
                MS2_ID=TestModel{i}.ImportantInf{J_Vec}.MS2_ID;
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_ID};
                [Max_val,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(MS2_ID,:));
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi01};
                [Max_val,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(MS2_ID,:));
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi02};
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[MS2_ID,Max_posi01,Max_posi02];
                
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};
            case 2 %%[0 1 0]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});  
                MS2_ID=TestModel{i}.ImportantInf{J_Vec}.MS2_ID;
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_ID};
                [Max_val,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(:,MS2_ID));
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi01};
                [Max_val,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(MS2_ID,:));
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi02};
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[Max_posi01,MS2_ID,Max_posi02];
                peptide{i} = Information_Matrix_Total{i, 9};
                protein{i} = Information_Matrix_Total{i, 14};
                iso{i} = Information_Matrix_Total{i, 16};
            case 4 %%[0 0 1]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID=TestModel{i}.ImportantInf{J_Vec}.MS2_ID;
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_ID};
                [Max_val,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(:,MS2_ID));
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi01};
                [Max_val,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(:,MS2_ID));
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi02};
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[Max_posi01,Max_posi02,MS2_ID];
                peptide{i} = Information_Matrix_Total{i, 17};
                protein{i} = Information_Matrix_Total{i, 22};
                iso{i} = Information_Matrix_Total{i, 24};
            case 3 %%[1 1 0]
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID01=TestModel{i}.ImportantInf{1}.MS2_ID;
                MS2_ID02=TestModel{i}.ImportantInf{2}.MS2_ID;
                
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_ID01};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(MS2_ID01,:));
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_ID02};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(MS2_ID02,:));
                
                if Max_val01>=Max_val02                
                    FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi01};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,MS2_ID02,Max_posi01];
                else                
                    FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{Max_posi02};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,MS2_ID02,Max_posi02];
                end
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};
                
            case 5 %%[1 0 1]
                
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID01=TestModel{i}.ImportantInf{1}.MS2_ID;
                MS2_ID02=TestModel{i}.ImportantInf{3}.MS2_ID;
                
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_ID01};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(MS2_ID01,:));
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_ID02};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(:,MS2_ID02));
                
                if Max_val01>=Max_val02                
                    FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi01};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,Max_posi01,MS2_ID02];
                else                
                    FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{Max_posi02};
                    FinalResult{i}.Cal_Posi=[MS2_ID01,Max_posi02,MS2_ID02];
                end
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};
                
            case 6 %%[0 1 1]
                
                [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel{i});
                MS2_ID01=TestModel{i}.ImportantInf{2}.MS2_ID;
                MS2_ID02=TestModel{i}.ImportantInf{3}.MS2_ID;
                
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_ID01};
                [Max_val01,Max_posi01]=max(ScoreMatrix01{ChooseScoreMatrix}(:,MS2_ID01));
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_ID02};
                [Max_val02,Max_posi02]=max(ScoreMatrix02{ChooseScoreMatrix}(:,MS2_ID02));
                
                if Max_val01>=Max_val02                
                    FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi01};
                    FinalResult{i}.Cal_Posi=[Max_posi01,MS2_ID01,MS2_ID02];
                else                
                    FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{Max_posi02};
                    FinalResult{i}.Cal_Posi=[Max_posi02,MS2_ID01,MS2_ID02];
                end
                FinalResult{i}.ScoreMatrix01=ScoreMatrix01;
                FinalResult{i}.ScoreMatrix02=ScoreMatrix02;
                FinalResult{i}.J_Vec=J_Vec;
                
                peptide{i} = Information_Matrix_Total{i, 9};
                protein{i} = Information_Matrix_Total{i, 14};
                iso{i} = Information_Matrix_Total{i, 16};
                
            case 7 %%[1 1 1]
                MS2_IDA=TestModel{i}.ImportantInf{1}.MS2_ID;
                MS2_IDB=TestModel{i}.ImportantInf{2}.MS2_ID;
                MS2_IDC=TestModel{i}.ImportantInf{3}.MS2_ID;                
                FinalResult{i}.ElutionProfile{1}=TestModel{i}.ImportantInf{1}.IntervalElutionMatrix{MS2_IDA};
                FinalResult{i}.ElutionProfile{2}=TestModel{i}.ImportantInf{2}.IntervalElutionMatrix{MS2_IDB};
                FinalResult{i}.ElutionProfile{3}=TestModel{i}.ImportantInf{3}.IntervalElutionMatrix{MS2_IDC};
                FinalResult{i}.J_Vec=J_Vec;
                FinalResult{i}.Cal_Posi=[MS2_IDA,MS2_IDB,MS2_IDC];      
                peptide{i} = Information_Matrix_Total{i, 1};
                protein{i} = Information_Matrix_Total{i, 6};
                iso{i} = Information_Matrix_Total{i, 8};
                
        end

end





num=1;
for i=1:length(FinalResult)    
    ID_nd=find(ID_nondetected-i==0);    
    if isempty(ID_nd)
        for j=1:3    
            IntensitySum=sum(FinalResult{i}.ElutionProfile{j},1);
            if length(IntensitySum)==8
                HLR(i,j)=sum(sum(IntensitySum(5:6)))/sum(sum(IntensitySum(1:2)));  
            else                
                HLR(i,j)=0;
            end
        end
    else        
        HLR(i,j)=0;
    end
end

HLR_Vector=HLR(:,1);
HLR_Vectorv1=HLR_Vector(HLR_Vector>=0 & HLR_Vector<=1000);
figure
hist(HLR_Vectorv1,[0:0.1:20])

nu=0;
ID_withKLabeled=[];
for i=1:size(Information_Matrix_Total,1)
    MZList01=Information_Matrix_Total{i,7};
    MZList02=Information_Matrix_Total{i,15};
    MZList03=Information_Matrix_Total{i,23};

    if length(MZList01)==8
        if MZList01(1)-MZList01(5)~=0
            nu=nu+1;ID_withKLabeled=[ID_withKLabeled;i];
        end
    else
            if length(MZList02)==8
                if MZList02(1)-MZList02(5)~=0
                    nu=nu+1;ID_withKLabeled=[ID_withKLabeled;i];
                end
            else
                    if length(MZList03)==8
                        if MZList03(1)-MZList03(5)~=0
                            nu=nu+1;ID_withKLabeled=[ID_withKLabeled;i];
                        end
                    end 
            end    
    end
end
save salicInfo HLR peptide protein FinalResult TestModel Information_Matrix_Total iso
load salicInfo
HLR_Vectorv2=HLR_Vector(ID_withKLabeled);


%test bench






 
