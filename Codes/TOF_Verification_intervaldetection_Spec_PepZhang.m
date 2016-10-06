function IntervalList=TOF_Verification_intervaldetection_Spec_PepZhang(iso,xics)

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
%maxGap=5;
mininterval=5;
klthreshold=-4;%% orbit -2.5 TOF -3.5
%Times_noise_std=3;%6;
%instrument_h_int=10^7;
%de_range=10^5;%10^4;
N_filter=5;

if iso(1)>iso(2)
   xiclist=[1 5];
else
   xiclist=[2 6];
end    

Weight=ones(N_filter,1)./N_filter;
[totalscans b]=size(xics);
XIC_afterfilter_1sttime=zeros(totalscans,b);
XIC_afterfilter_2ndtime=zeros(totalscans,b);
XIC_afterfilter_3ndtime=zeros(totalscans,b);
for j=1:6   
    XIC_afterfilter_1sttime(1:totalscans,j)=filter2(Weight,xics(:,j)); 
    XIC_afterfilter_2ndtime(1:totalscans,j)=filter2(Weight,XIC_afterfilter_1sttime(:,j)); 
    XIC_afterfilter_3ndtime(1:totalscans,j)=filter2(Weight,XIC_afterfilter_2ndtime(:,j)); 
end    

XICforIntervalDetection=xics(:,xiclist(1))+xics(:,xiclist(2));
SmoothXICforIntervalDetection=filter2(Weight, XICforIntervalDetection); 
SmoothXICforIntervalDetection=filter2(Weight, SmoothXICforIntervalDetection);
SmoothXICforIntervalDetection=filter2(Weight, SmoothXICforIntervalDetection);
th=getNoiseThreshold(SmoothXICforIntervalDetection,XICforIntervalDetection,noiseThresholdLevel);
diffXIC=SmoothXICforIntervalDetection;
diffXIC(2:end)=SmoothXICforIntervalDetection(2:end)-SmoothXICforIntervalDetection(1:end-1);
[intervalList, intervalCount]=getInterval(SmoothXICforIntervalDetection,th,mininterval);

IntervalList.intervallistAfterInitialDetection=intervalList;
[intervalList,intervalCount]=splitInterval(SmoothXICforIntervalDetection,diffXIC,intervalList,intervalCount, mininterval);
[IntervalList.intervallistAfterSplit totalInterval]=removeShortIntervals(intervalList,mininterval);
% ms2_time=Orbi_ms2information;


%%%%%%%%%%%%%%% interval boundary detection and   long's
%%%%%%%%%%%%%%% criteria
Total_interval=IntervalList.intervallistAfterSplit;
Pep_iso=iso;
if totalInterval==0
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
       IntervalList.intervallist_afterBoundaryCheck(intervalIndex,:)=[Scan_start_afterfilter,Scan_end_afterfilter];  
        % check if the length is less than four, set the corresponding interval
        %indicator as zero
       if Scan_end_afterfilter-Scan_start_afterfilter+1<=mininterval
           intervalIndicator(intervalIndex)=0;
       end
       
       if  intervalIndicator(intervalIndex)==1;
         [Judge, Labelling_efficiency]=checkDataModelO18(Pep_iso,intervalsdata01v2,klthreshold);
          if Judge==1
             Judge=checkPeakShapes(intervalsdata01v2,2);
          end     
                   
          if Judge==0
             IntervalList.intervallist_afterBoundaryCheck(intervalIndex,:)=[Scan_start_afterfilter,Scan_start_afterfilter];                   
           end 
       end 

    end    
    
end    
if totalInterval==0
     IntervalList.intervallist_after7_combine=[0 0];
     IntervalList.totalInterval=0;
else     
    
 [Total_intervallist_new totalInterval]=removeShortIntervals(IntervalList.intervallist_afterBoundaryCheck,mininterval);
 
 if totalInterval>0
     IntervalList.intervallist_after7_combine=Total_intervallist_new;
     IntervalList.totalInterval=totalInterval;    
 else
      IntervalList.intervallist_after7_combine=[0 0];
     IntervalList.totalInterval=0;
 end    
end

%%%%%%%%%%%%%%%


%   colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
% % 
%      figure; hold on;
%      for j=1:6
%          plot(xics(:,j),colorarray(j))
% %         Threshold(j)=getThreshold_range(XICs_Orbi01{j}(:,ID),10000000, 10000,Times_noise_std);
% %         hold on
% %         plot(1:length(XICs_Orbi01{j}(:,ID)),Threshold(j)*ones(1,length(XICs_Orbi01{j}(:,ID))),colorarray(j))
%      end
% %     
%      height=5*th;
%     stem(IntervalList.intervallistAfterInitialDetection(:,1),height*ones(length(IntervalList.intervallistAfterInitialDetection(:,1)),1),'ro')
%     stem(IntervalList.intervallistAfterInitialDetection(:,2),height*ones(length(IntervalList.intervallistAfterInitialDetection(:,1)),1),'ko')
%     stem(IntervalList.intervallistAfterSplit(:,1),1.5*height*ones(length(IntervalList.intervallistAfterSplit(:,1)),1),'r+')
%     stem(IntervalList.intervallistAfterSplit(:,2),1.5*height*ones(length(IntervalList.intervallistAfterSplit(:,1)),1),'k+')
%     stem(IntervalList.intervallist_afterBoundaryCheck(:,1),2*height*ones(length(IntervalList.intervallist_afterBoundaryCheck(:,1)),1),'r>')
%     stem(IntervalList.intervallist_afterBoundaryCheck(:,2),2*height*ones(length(IntervalList.intervallist_afterBoundaryCheck(:,1)),1),'k>')
%     ll=length(IntervalList.intervallist_afterBoundaryCheck(:,1));
%     stem(IntervalList.intervallist_after7_combine(:,1),3*height*ones(length(IntervalList.intervallist_after7_combine(:,1)),1),'rs')
%     stem(IntervalList.intervallist_after7_combine(:,2),3*height*ones(length(IntervalList.intervallist_after7_combine(:,1)),1),'ks')
%    
%     for tt=1:length(IntervalList.intervallist_after7_combine(:,1))
%         text(IntervalList.intervallist_after7_combine(tt,1),3*height,num2str(tt)); 
%          text(IntervalList.intervallist_after7_combine(tt,2),3*height,num2str(tt)); 
%     end    
%    
    

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

