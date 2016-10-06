%%%%%%%%%%%%%%%%%% read Maxquant txt file and generate Information matrix
%%%%%%%%%%%%%%% Total peptide information
%%%%%%%%%%%%%%% column name: peptide_seq cs mass  ms2scannumber HLR protein
%%%%%%%%%%%%%%% totalmzlist iso (1-8 column for 1st replicate 9-16 for 2nd 17-24 for 3rd)
%%%%%% note: load Maxquant processing result of Pre-ratio replicates (3 data with same HLR)
%%%%%% the sum of iso is not normalized to 1.

if runPreprocessing~=0
    
    disp('Reading MS2 information from maxquant......')
    [Information_Matrix_K_labeled_A, Information_Matrix_A]=readMaxquantResult(path1, PEPTIDE, MSMS, OXIDATION, LABEL_TECH);
    [Information_Matrix_K_labeled_B, Information_Matrix_B]=readMaxquantResult(path2, PEPTIDE, MSMS, OXIDATION, LABEL_TECH);
    [Information_Matrix_K_labeled_C, Information_Matrix_C]=readMaxquantResult(path3, PEPTIDE, MSMS, OXIDATION, LABEL_TECH);
    
    %%%%%%%%%% combine information matrix from 3 input matrix (Remove redundant entries)
    disp('MS2 information Combining......')
    Information_Matrix_Total=CombineOriginalInformationMatrix(Information_Matrix_K_labeled_A,Information_Matrix_K_labeled_B);
    Information_Matrix_Total=CombineOriginalInformationMatrix(Information_Matrix_Total,Information_Matrix_K_labeled_C);
    
    [nUnion, ~] = size(Information_Matrix_Total);
    [nData1, ~] = size(Information_Matrix_K_labeled_A);
    [nData2, ~] = size(Information_Matrix_K_labeled_B);
    [nData3, ~] = size(Information_Matrix_K_labeled_C);
    %%%%%%%%%%%%%%%%%%
    ms2IdentifiedId = [cell2mat(Information_Matrix_Total(:,2)), cell2mat(Information_Matrix_Total(:,10)), cell2mat(Information_Matrix_Total(:,18))];
    allIdentifiedId = find(ms2IdentifiedId(:,1)>0&ms2IdentifiedId(:,2)>0|ms2IdentifiedId(:,1)>0&ms2IdentifiedId(:,3)>0);
    Information_Matrix_Total = Information_Matrix_Total(allIdentifiedId, :);
    
%     Information_Matrix_Total = Information_Matrix_Total(1:500,:);
    %%%%%%%%%%%%%%%%% read mzxml data
    

    [retentiontl1_A,~,peakl_A,retentiont_A,MZInt_l1l2_A]=readrawdata(MzxmlFile_path1);

    
    [retentiontl1_B,~,peakl_B,retentiont_B,MZInt_l1l2_B]=readrawdata(MzxmlFile_path2);

    
    [retentiontl1_C,~,peakl_C,retentiont_C,MZInt_l1l2_C]=readrawdata(MzxmlFile_path3);

    
    
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
    toc
    XICs_B=getXIC_LC_new(peakl_B,TOF_totalmzList,tolerance);
    XICs_C=getXIC_LC_new(peakl_C,TOF_totalmzList,tolerance);
    clear peakl_A peakl_B peakl_C

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% interval detection
    %%%%%%%%%%%%%%%%% the returned data structure contains intervals of LC peaks. the
    %%%%%%%%%%%%%%%%% final intervallist is called intervallist_after7_combine
    
    pep_number=size(XICs_A,2)/8;
    
    IntervalListBig = cell(pep_number, 3);
    for i=1:pep_number
        disp(['detecting ', num2str(i), 'th peptide. total: ', num2str(pep_number), '.'])
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
            [IntervalList,Jud_detect_good]=SILAC_Verification_intervaldetection_Spec_Pepv1(isoList,pepXICs, kl_thresHold);
            IntervalListBig{i,filecount}=IntervalList;
            IntervalListBig{i,filecount}.chargeState=cs;
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
        INTERVAL_Matrix(i, :)
    end
    
    % diagnostic code

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
    save([save_path,ExperimentID, '_',  'XICandIntervalInformation2'], 'XICs_A', 'XICs_B', 'XICs_C', 'IntervalListBig', 'INTERVAL_Matrix', 'Information_Matrix_Total', 'nData1', 'nData2', 'nData3', 'nUnion')
    save([save_path,ExperimentID, '_', 'RetentionTime'], 'retentiontl1_A', 'retentiontl1_B', 'retentiontl1_C', 'retentiont_A', 'retentiont_B', 'retentiont_C')

else
    load([save_path,ExperimentID, '_',  'XICandIntervalInformation2'], 'XICs_A', 'XICs_B', 'XICs_C', 'IntervalListBig', 'INTERVAL_Matrix', 'Information_Matrix_Total', 'nData1', 'nData2', 'nData3', 'nUnion')
    load([save_path,ExperimentID, '_', 'RetentionTime'], 'retentiontl1_A', 'retentiontl1_B', 'retentiontl1_C', 'retentiont_A', 'retentiont_B', 'retentiont_C')
end
    
    
