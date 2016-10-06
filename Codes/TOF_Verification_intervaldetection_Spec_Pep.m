function IntervalList=TOF_Verification_intervaldetection_Spec_Pep(iso,xics)

%%%%%%%%%% this function is detect interval for one peptide and judge which
%%%%%%%%%% interval is ok for it after 1~5steps and KL interval scan
%%%%%%%%%% detection and long's criteria and combine with O16 and O18
%%%%%%%%%% intervals


% Pepcommon0102v1=pep01;
% iso=iso01;
% Orbi_ms2information01v1=ms2information01;
% XICs_Orbi01=XICs_Orbi01;
% Orbi_retentiont01l1=retentiont01l1;

%%%%%%%%%%%%% Original interval detection
%%%%%%%%%%%%% 1. Auto threshold set up
%%%%%%%%%%%%% 2. Smooth the XIC before the detection
%%%%%%%%%%%%% 3. Find the point with max intensity
%%%%%%%%%%%%% 4. Get the first order derivatives
%%%%%%%%%%%%% 5. Find the two smallest derivatives points around
noiseThresholdLevel=3;
maxGap=5;
mininterval=5;
klthreshold=-2.5;
Times_noise_std=3;%6;
instrument_h_int=10^7;
de_range=10^5;%10^4;
N_filter=5;
Weight=ones(N_filter,1)./N_filter;
[totalscans b]=size(xics);
XIC_afterfilter_1sttime=zeros(totalscans,b);
XIC_afterfilter_2ndtime=zeros(totalscans,b);
XIC_afterfilter_3ndtime=zeros(totalscans,b);
for j=1:6   
    XIC_afterfilter_1sttime(1:totalscans,j)=filter2(Weight,xics(:,j)); 
    XIC_afterfilter_2ndtime(1:totalscans,j)=filter2(Weight,XIC_afterfilter_1sttime(:,j)); 
    XIC_afterfilter_3ndtime(1:totalscans,j)=filter2(Weight,XIC_afterfilter_2ndtime(:,j)); 
    th(j)=getNoiseThreshold(XIC_afterfilter_3ndtime(1:totalscans,j),xics(:,j),noiseThresholdLevel);
end    
TotalTH=sum(th);
% ms2_time=Orbi_ms2information;

if iso(1)>iso(2)
   xiclist=[1 5];
else
   xiclist=[2 6];
end    
diffXIC=zeros(1,totalscans);
for xiccount=1:2
    j=xiclist(xiccount);
    [intervalList, intervalCount]=getInterval(XIC_afterfilter_3ndtime(:,j),th(j),mininterval);
   % intervallist{xiccount}=intervalList;
    diffXIC(1)=XIC_afterfilter_3ndtime(1,j);
    diffXIC(2:end)=XIC_afterfilter_3ndtime(2:end,j)-XIC_afterfilter_3ndtime(1:end-1,j);
    [intervalList,intervalCount]=splitInterval(XIC_afterfilter_3ndtime(:,j),diffXIC,intervalList,intervalCount, mininterval);
     intervallist_after5{xiccount}=intervalList;         
     intervalCount_vector(xiccount)=intervalCount;      
 end

%IntervalList.intervallist=intervallist;%{Record(:,1)}(Record(P_m,2),:);
IntervalList.intervallist_after5=intervallist_after5;
Total_intervallist=[intervallist_after5{1}; intervallist_after5{2}];
[IntervalList.intervallistv1 totalInterval]=removeShortIntervals(Total_intervallist,mininterval);




%%%%%%%%%%%%%%% interval boundary detection and   long's
%%%%%%%%%%%%%%% criteria
Total_interval=IntervalList.intervallistv1;
Pep_iso=iso;
if totalInterval==0
    IntervalList.intervallist_after7=[0 0];
    IntervalList.Labelling_effeciency_after7= [0 0];
     IntevalList.totalInterval=0;
    
else
    intervalIndicator=ones(totalInterval,1);
    Labelling_effeciency=zeros(totalInterval,1);
    for intervalIndex=1:totalInterval
        Scan_start=Total_interval(intervalIndex,1);
        Scan_end=Total_interval(intervalIndex,2);
 
       intervalsdata01v1_afterfilter=XIC_afterfilter_3ndtime(Scan_start:Scan_end,:);
       LC_profile=sum(intervalsdata01v1_afterfilter,2);
       [temp, Ms2_scannumber]=max(LC_profile);
       Ms2_scannumber=Ms2_scannumber+Scan_start-1;

       
       [intervalsdata01v2 ,Scan_start_afterfilter,Scan_end_afterfilter]=...
        BoundaryDetection(intervalsdata01v1_afterfilter, Scan_start,Ms2_scannumber,klthreshold);  
       IntervalList.intervallist_after6(intervalIndex,:)=[Scan_start_afterfilter,Scan_end_afterfilter];  
        % check if the length is less than four, set the corresponding interval
        %indicator as zero
       if Scan_end_afterfilter-Scan_start_afterfilter+1<=mininterval
           intervalIndicator(intervalIndex)=0;
       end
       
       if  intervalIndicator(intervalIndex)==1;
         [Judge, Labelling_efficiency]=checkDataModelO18(Pep_iso,intervalsdata01v2,klthreshold);
        
         intervalIndicator(intervalIndex)=Judge;
       end 

    end    


    ID_keep=(intervalIndicator==1);
    if isempty(ID_keep)==0
       IntervalList.intervallist_after7= IntervalList.intervallist_after6(ID_keep,:);
    
        IntevalList.totalInterval=sum(ID_keep);
           
    else   
        IntervalList.intervallist_after7=[0 0];

         IntevalList.totalInterval=0;
        
    end 
end    




%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% combine O16 and O18 intervals based on iso mono and 1st
if IntevalList.totalInterval>0
        Total_intervallist=IntervalList.intervallist_after7;
  
        %%%%%%%%%%%%% combine the O16 and O18 intervals
      
       [Total_intervallist_sort,totalInterval]=clearOverlappingPeaks(Total_intervallist);        
        [Total_intervallist_new totalInterval]=mergeClosePeaks(Total_intervallist_sort, XIC_afterfilter_3ndtime, klthreshold-2, maxGap);
        [Total_intervallist_new totalInterval]=removeShortIntervals(Total_intervallist_new,mininterval);
        if totalInterval>0
            IntervalList.intervallist_after7_combine=Total_intervallist_new;
             IntervalList.totalInterval=totalInterval;
        else 
            IntervalList.intervallist_after7_combine=[0 0];
             IntervalList.totalInterval=0;
        end    
else       
       
    IntervalList.intervallist_after7_combine=[0 0];  
    IntervalList.intervallist_after6_combine=[0 0]; 
     IntervalList.totalInterval=0;
end            

% Re-adjust the boundary after combining information
if  IntervalList.totalInterval>0
    for k=1: IntervalList.totalInterval
        Scan_data=IntervalList.intervallist_after7_combine(k,:);
        if Scan_data(1)~=0 && Scan_data(2)~=0
               intervalsdata_afterfilter=XIC_afterfilter_3ndtime(Scan_data(1):Scan_data(2),:);
               LC_profile=sum(intervalsdata_afterfilter,2);
               [MaxVal MaxScan]=max(LC_profile);
               MaxScan=MaxScan+Scan_data(1)-1;
                [intervalsdata ,Scan_start_afterfilter,Scan_end_afterfilter]=BoundaryDetection(intervalsdata_afterfilter, Scan_data(1),MaxScan,klthreshold);  
                IntervalList.intervallist_after7_combine(k,:)=[Scan_start_afterfilter,Scan_end_afterfilter];  
                 if Scan_end_afterfilter-Scan_start_afterfilter+1>=mininterval
                   [Judge, Labelling_efficiency]=checkDataModelO18(Pep_iso,intervalsdata,klthreshold-1);
                    if Judge==1
                       Judge=checkPeakShapes(intervalsdata,0.5);
                    end     
                   
                   if Judge==0
                         IntervalList.intervallist_after7_combine(k,:)=[Scan_start_afterfilter,Scan_start_afterfilter];                   
                    end 
                 end   
        end
    end
     [Total_intervallist_new totalInterval]=removeShortIntervals(IntervalList.intervallist_after7_combine,mininterval);
     if totalInterval>0
         IntervalList.intervallist_after7_combine=Total_intervallist_new;
         IntervalList.totalInterval=totalInterval;
     else
          IntervalList.intervallist_after7_combine=[0 0];
         IntervalList.totalInterval=0;
     end    
         
end
%%%%%%%%%%%%%% Recalibrate for the boundrys calculate the labelling effeciency for each final intervals

%colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
if   IntervalList.totalInterval>0
 for k=1:  IntervalList.totalInterval
    Scan_data=IntervalList.intervallist_after7_combine(k,:);
    if Scan_data(1)~=0 && Scan_data(2)~=0
        
           intervalsdata=xics(Scan_data(1):Scan_data(2),:);
           [AbundanceO16O18 Labelling_efficiency]=estimateO16O18ModifiedYao(iso,intervalsdata);
            IntervalList.Labelling_efficiency_after7_combine(k)=Labelling_efficiency;
      %  intensity=sum(Interval_data,1);
    %    est=O18rateLinear(iso(1:6),intensity,minf,maxf,rangerate,1);
     %   IntervalList.Labelling_efficiency_after7_combine(k)=est(3);
        
    else
        IntervalList.Labelling_efficiency_after7_combine(k)=0;        
    end
    
  end
end%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% check
% ID_only_ginterval=[];
% ID_long_Judge=[];
% for i=1:length(ID_Good_groundinterval_correction)
%     ID=ID_Good_groundinterval_correction(i);    
%     IntervalList{ID}.intervallist{1};
%     IntervalList{ID}.intervallistv1;
%     ID_good_judge=find(IntervalList{ID}.Judge==1);
%     if ~isempty(ID_good_judge)
%         ID_long_Judge=[ID_long_Judge;ID, length(ID_good_judge)];
%     end
%     if size(IntervalList{ID}.intervallistv1,1)==1
%         ID_only_ginterval=[ID_only_ginterval; i, ID];        
%     end   
% end
% 
% length(ID_Good_groundinterval_correction)
% length(ID_only_ginterval)
% length(ID_long_Judge)

% for i=100:105%:210;
%     ID=ID_Good_groundinterval_correction(i);
%     colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
% 
%     figure
%     for j=1:6
%         plot(XICs_Orbi01{j}(:,ID),colorarray(j))
%         Threshold(j)=getThreshold_range(XICs_Orbi01{j}(:,ID),10000000, 10000,Times_noise_std);
%         hold on
%         plot(1:length(XICs_Orbi01{j}(:,ID)),Threshold(j)*ones(1,length(XICs_Orbi01{j}(:,ID))),colorarray(j))
%     end
%     
%     height=5*max(Threshold);
%     stem(IntervalList{ID}.intervallist{1}(:,1),height*ones(length(IntervalList{ID}.intervallist{1}(:,1)),1),'ro')
%     stem(IntervalList{ID}.intervallist{1}(:,2),height*ones(length(IntervalList{ID}.intervallist{1}(:,1)),1),'ko')
%     stem(IntervalList{ID}.intervallist_after5{1}(:,1),1.5*height*ones(length(IntervalList{ID}.intervallist_after5{1}(:,1)),1),'r*')
%     stem(IntervalList{ID}.intervallist_after5{1}(:,2),1.5*height*ones(length(IntervalList{ID}.intervallist_after5{1}(:,1)),1),'k*')
%     stem(IntervalList{ID}.intervallist_after6(:,1),2*height*ones(length(IntervalList{ID}.intervallist_after6(:,1)),1),'rs')
%     stem(IntervalList{ID}.intervallist_after6(:,2),2*height*ones(length(IntervalList{ID}.intervallist_after6(:,1)),1),'ks')
%     stem(IntervalList{ID}.intervallist_after7(:,1),2.5*height*ones(length(IntervalList{ID}.intervallist_after7(:,1)),1),'r+')
%     stem(IntervalList{ID}.intervallist_after7(:,2),2.5*height*ones(length(IntervalList{ID}.intervallist_after7(:,1)),1),'k+')
%     stem(IntervalList{ID}.intervallist_after7_combine(:,1),3*height*ones(length(IntervalList{ID}.intervallist_after7_combine(:,1)),1),'rd')
%     stem(IntervalList{ID}.intervallist_after7_combine(:,2),3*height*ones(length(IntervalList{ID}.intervallist_after7_combine(:,1)),1),'kd')
% 
%     
% %     ID_long=find(IntervalList{ID}.Judge==1);
% %     stem(IntervalList{ID}.intervallist_after6(ID_long,1),2.5*10^7*ones(length(IntervalList{ID}.intervallistv1(ID_long,1)),1),'r+')
% %     stem(IntervalList{ID}.intervallist_after6(ID_long,2),2.5*10^7*ones(length(IntervalList{ID}.intervallistv1(ID_long,2)),1),'k+')
%     grid on
% 
% end

% ID_after7_Judgenotinterval=[];
% for j=1:length(ID_Good_groundinterval_correction)
%     ID=ID_Good_groundinterval_correction(j);
%     if sum(sum(IntervalList{ID}.intervallist_after7_combine))==0
%         ID_after7_Judgenotinterval=[ID_after7_Judgenotinterval;j ID];
%     end 
% end

% IntervalList{ID}.intervallist{5}
% IntervalList{ID}.intervallist_after5
% IntervalList{ID}.intervallist_after6
% ID_long=find(IntervalList{ID}.Judge==1);


% %%%%%%%%%%%%% QTOF same peptide XIC
% filename=['D:\Program\QTOF_replicate_identification\ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_01_f4.mzXML.centroid1.peak.mat'];
% load(filename);
% peakl01=peakl; 
% TOF_retentiont01l1=retentiont;
% clear peakl filename retentiont;
% 
% filename=['D:\Program\QTOF_replicate_identification\ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_02_f4.mzXML.centroid1.peak.mat'];
% load(filename);
% peakl02=peakl;
% TOF_retentiont02l1=retentiont;
% clear peakl filename retentiont;
% 
% filename=['D:\Program\QTOF_replicate_identification\ZHA_27_5574_21APR11_CELL_VEL_HUM_TT_JL_2D_03_f4.mzXML.centroid1.peak.mat'];
% load(filename);
% peakl03=peakl;
% TOF_retentiont03l1=retentiont;
% clear peakl filename retentiont;
% 
% for kkkk=100:105
%  ID=ID_Good_groundinterval_correction(kkkk);
% 
% TOF_spec_Pep=pep01{ID};
% TOF_spec_Pep_cs=chargestate01(ID);
% TOF_spec_Pep_mass=mass01(ID);
% TOF_spec_Pep_iso=iso01(ID,:);
% TOF_spec_Pep_mzlist=totalmzList01(ID,:);
% tolerance=20;
% for i=1:6
%     XICs_TOF01{i}=getXICs(peakl01,TOF_spec_Pep_mzlist(i),tolerance);
% end
% 
% mininterval=4;
% maxnointervals=20;
% 
% for i=1:length(TOF_spec_Pep_cs)
%      
%         xic=zeros(size(XICs_TOF01{1},1),6);
%         for j=1:6
%             xic(:,j)=XICs_TOF01{j}(:,i);
%             %%%%%%%interval detection: 
%             %%%%%%%1. threshold up and length of interval
%             %%%%%%%larger than 5 points
%             min_interval_length=6;
%             %%%%%%%%%%%% mean filter the XIC
%             N_filter=5;
%             Weight=ones(N_filter,1)./N_filter;
%             XIC_afterfilter_1sttime=filter2(Weight,xic(:,j)); 
%             XIC_afterfilter_2ndtime=filter2(Weight,XIC_afterfilter_1sttime); 
% %             figure
% %             plot(xic(:,j),'r')
% %             hold on
% %             plot(XIC_afterfilter_1sttime,'k')
% %             plot(XIC_afterfilter_2ndtime,'g')
%             %%%%%%%%%%%%
%             TOF_intervallist{j}=intervaldetectionv1(xic(:,j),XIC_afterfilter_2ndtime,mininterval,maxnointervals,min_interval_length);
%             %%%%%%%
%             if TOF_intervallist{j}(1,1)~=0 && TOF_intervallist{j}(1,2)~=0
%                 for k=1:size(TOF_intervallist{j},1)
%                     Scan_start01=TOF_intervallist{j}(k,1);
%                     Scan_end01=TOF_intervallist{j}(k,2);
%                     start01=TOF_retentiont01l1(Scan_start01);
%                     end01=TOF_retentiont01l1(Scan_end01);
%                     %%%%%2. smooth all the intervals
%                     N_filter_vector=1:2:2000;
%                     Interval_XIC=xic(Scan_start01:Scan_end01,j);
%                     P_filter=find(N_filter_vector<=round(length(Interval_XIC)/3));
%                     N_filter=N_filter_vector(P_filter(end));
%                     Weight=ones(N_filter,1)./N_filter;
%                     Interval_XIC_afterfilter_1sttime=filter2(Weight,Interval_XIC); 
%                     Interval_XIC_afterfilter=filter2(Weight,Interval_XIC_afterfilter_1sttime); 
%                     
% %                     figure
% %                     plot(Interval_XIC,'r');hold on;plot(Interval_XIC_afterfilter,'k')
%                     %%%%%
%                     %%%%%3. find the max peak point
%                     %%%%%4. find order one derivitative
%                     %%%%%5. find two smallest derivitative points around the max point
%                     [max_value, peak_start, peak_end]=Get_one_peak(Interval_XIC_afterfilter);
%                     %%%%%
%                     TOF_intervallist_after5{j}(k,:)=[Scan_start01+peak_start-1, Scan_start01+peak_end-1];
%                     %%%%%
%  
%                 end
%             end
%         end
%            
% end
% 
% %%%%%%%%%%%% combine O16 and O18 intervals based on iso mono and 1st
% for j=1:length(TOF_spec_Pep_cs)
%     
%     if TOF_spec_Pep_iso(j,1)>=TOF_spec_Pep_iso(j,2)
%             %%%%%%%%%%%%% combine the O16 and O18 intervals
%             ID_iso1=1;
%             ID_iso2=5;
%             pepXICs=zeros(size(XICs_TOF01{1},1),6);
%             for i=1:6
%                 pepXICs(:,i)=XICs_TOF01{i}(:,j);
%             end          
%             TOF_intervallistv1=[[TOF_intervallist{ID_iso1},ID_iso1*ones(size(TOF_intervallist{ID_iso1},1),1)];[TOF_intervallist{ID_iso2},ID_iso2*ones(size(TOF_intervallist{ID_iso2},1),1)]];
%             Total_intervallist_new=combine_O16_O18_intervals(TOF_intervallistv1,ID_iso1,ID_iso2,pepXICs);                
% %             Total_intervallist_new;
%     else
%             %%%%%%%%%%%%% combine the O16 and O18 intervals
%             ID_iso1=2;
%             ID_iso2=6;
%             pepXICs=zeros(size(XICs_TOF01{1},1),6);
%             for i=1:6
%                 pepXICs(:,i)=XICs_TOF01{i}(:,j);
%             end            
%             TOF_intervallistv1=[[TOF_intervallist{ID_iso1},ID_iso1*ones(size(TOF_intervallist{ID_iso1},1),1)];[TOF_intervallist{ID_iso2},ID_iso2*ones(size(TOF_intervallist{ID_iso2},1),1)]];
%             Total_intervallist_new=combine_O16_O18_intervals(TOF_intervallistv1,ID_iso1,ID_iso2,pepXICs); 
%     end
%    
% end            
% 
% ID_zeros=sum(Total_intervallist_new(:,1:2),2)==0;
% Total_intervallist_new(ID_zeros,:)=[];
% 
% %%%%%%%%%%%%%%% interval spearman correlation correction and long's
% %%%%%%%%%%%%%%% criteria
% for i=1:length(TOF_spec_Pep_cs)
%     TOF_Total_interval=Total_intervallist_new;
%     TOF_Pep_iso=TOF_spec_Pep_iso;
%     TOF_instrument_h_int=10^6;
%     TOF_de_range=10^4;
%     for j=1:6
%         TOF_XIC(:,j)=XICs_Orbi01{j}(:,i);
%         TOF_Interval_int_threshold(j)=getThreshold_range(TOF_XIC,TOF_instrument_h_int, TOF_de_range,3);
%     end
%     
%     for k=1:size(TOF_Total_interval,1)        
%         Scan_start=TOF_Total_interval(k,1);
%         Scan_end=TOF_Total_interval(k,2);
%         ID_iso=TOF_Total_interval(k,3);
%         [V_max,P_max]=max(XICs_TOF01{ID_iso}(Scan_start:Scan_end,i));
%         TOF_intervaldata=zeros(Scan_end-Scan_start+1,6);
%         TOF_intervalsdata01v1_afterfilter=zeros(Scan_end-Scan_start+1,6);
%         N_mean_filter=5;
%         for j=1:6
%             TOF_intervaldata(:,j)=XICs_TOF01{j}(Scan_start:Scan_end,i);
%             MS2_Peak_Intensity01(j)=TOF_intervaldata(P_max,j);
%             
%             %%%%%%%%%%%%%%%%%% XIC processed by mean filter 
%             Weight=ones(N_mean_filter,1)./N_mean_filter;
%             TOF_intervalsdata01v1_afterfilter(:,j)=filter2(Weight,TOF_intervaldata(:,j)); 
%             MS2_Peak_Intensity01_afterfilter(j)=TOF_intervalsdata01v1_afterfilter(P_max,j);
%             %%%%%%%%%%%%%%%%%%                      
%         end
%         Ms2_scannumber=P_max+Scan_start-1;
%         [TOF_intervalsdata01v2 ,Scan_start_afterfilter,Scan_end_afterfilter]=MS2intervalfilter(TOF_intervaldata, MS2_Peak_Intensity01, TOF_intervalsdata01v1_afterfilter, MS2_Peak_Intensity01_afterfilter, Scan_start, Scan_end, Ms2_scannumber);
%         TOF_intervallist_after6(k,:)=[Scan_start_afterfilter,Scan_end_afterfilter];
%         %%%%%%%%%%%%%%%%%%%%%%%%%%long's peptide checks criteria
% 
%         Judge=long_criteria(TOF_Pep_iso,intervalsdata01v2,TOF_Interval_int_threshold);
% %         Scan_start_afterfilter;
% %         Scan_end_afterfilter;
%         %%%%%%%%%%%%%
%         AA(k)=Judge;
%     
%     end
%     
% end
% %%%%%%%%%%%%%%%
% 
% for i=1
%     colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
% 
%     figure
%     for j=1:6
%         plot(XICs_TOF01{j}(:,i),colorarray(j))
%         Threshold(j)=getThreshold_range(XICs_TOF01{j}(:,i),1000000, 100000, 3);
%         hold on
%         plot(1:length(XICs_TOF01{j}(:,i)),Threshold(j)*ones(1,length(XICs_TOF01{j}(:,i))),colorarray(j))
%     end
%     height=5*max(Threshold);
%     stem(Total_intervallist_new(:,1),height*ones(length(Total_intervallist_new(:,1)),1),'ro')
%     stem(Total_intervallist_new(:,2),height*ones(length(Total_intervallist_new(:,1)),1),'ko')
%     stem(TOF_intervallist_after5{1}(:,1),1.5*height*ones(length(TOF_intervallist_after5{1}(:,1)),1),'r*')
%     stem(TOF_intervallist_after5{1}(:,2),1.5*height*ones(length(TOF_intervallist_after5{1}(:,1)),1),'k*')
%     stem(TOF_intervallist_after6(:,1),1.5*height*ones(length(TOF_intervallist_after6(:,1)),1),'rs')
%     stem(TOF_intervallist_after6(:,2),1.5*height*ones(length(TOF_intervallist_after6(:,1)),1),'ks')
% %     ID_long=find(TOF_Judge==1);
% %     stem(TOF_intervallist_after6(ID_long,1),2.5*10^7*ones(length(TOF_intervallistv1(ID_long,1)),1),'r+')
% %     stem(TOF_intervallist_after6(ID_long,2),2.5*10^7*ones(length(TOF_intervallistv1(ID_long,2)),1),'k+')
%     grid on
% 
% end
% 
% 
% colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
% 
% % figure
% % for j=1:6
% %     plot(XICs_TOF01{j}(:,i),colorarray(j))
% %     Threshold(j)=getThreshold_range(XICs_TOF01{j}(:,i),10000, 2000,3);
% %     hold on
% %     plot(1:length(XICs_TOF01{j}(:,i)),Threshold(j)*ones(1,length(XICs_TOF01{j}(:,i))),colorarray(j))
% % end
% % grid on;
% 
% clear TOF_intervallist TOF_intervallist_after5 TOF_intervallist_after6
% 
% end
% %%%%%%%%%%%%%




