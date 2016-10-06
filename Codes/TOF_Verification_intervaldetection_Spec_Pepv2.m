function [IntervalList,Jud_detect_good]=TOF_Verification_intervaldetection_Spec_Pepv2(Pep,iso,XICs)

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
mininterval=4;
maxnointervals=20;
Times_noise_std=3;%6;
instrument_h_int=10^7;
de_range=10^5;%10^4;

% ms2_time=Orbi_ms2information;
xic=XICs;


for j=1:6
%     xic(:,j)=XICs(:,i);
    %%%%%%%interval detection: 
    %%%%%%% 1. threshold up and length of interval
    %%%%%%%     larger than 5 points
    min_interval_length=6;
    %%%%%%% 2.Smooth the XIC before the detection / mean filter the XICs
    N_filter=5;
    Weight=ones(N_filter,1)./N_filter;
    XIC_afterfilter_1sttime=filter2(Weight,xic(:,j)); 
    XIC_afterfilter_2ndtime=filter2(Weight,XIC_afterfilter_1sttime); 
%             figure
%             plot(xic(:,j),'r')
%             hold on
%             plot(XIC_afterfilter_1sttime,'k')
%             plot(XIC_afterfilter_2ndtime,'g')
%             grid on
%             plot(1:length(XIC_afterfilter_1sttime),th*ones(1,length(XIC_afterfilter_1sttime)))
%             Height=600;
%             stem((intervallist{j}(:,1)),2*Height*ones(length(intervallist{j}(:,1)),1),'rd')
%             stem((intervallist{j}(:,2)),2*Height*ones(length(intervallist{j}(:,2)),1),'kd')
    
    %%%%%%%%%%%%
    
    %%%%%%% interval detection based on the XICs after mean filter
    %%%%%%% smoothing

    % th=getThreshold(xic);
    th=TOF_getThreshold_range(xic(:,j),instrument_h_int, de_range,Times_noise_std);
    intervallist{j}=TOF_intervaldetectionv1(xic(:,j),XIC_afterfilter_2ndtime,mininterval,maxnointervals,min_interval_length,th);
    %%%%%%%
    
    %%%%%%%    
    if intervallist{j}(1,1)~=0 && intervallist{j}(1,2)~=0
        num=1;
        for k=1:size(intervallist{j},1)
            Scan_start01=intervallist{j}(k,1);
            Scan_end01=intervallist{j}(k,2);
%             start01=Orbi_retentiontl1(Scan_start01);
%             end01=Orbi_retentiontl1(Scan_end01);
%             if ms2_time>=start01 && ms2_time<=end01
%                 Record=[Record;j,k];
%                 MaxValue=[MaxValue;max(xic(Scan_start01:Scan_end01,j))];
%             end
            %%%%%2. smooth all the intervals
            N_filter_vector=1:2:2000;
            Interval_XIC=xic(Scan_start01:Scan_end01,j);
            P_filter=find(N_filter_vector<=round(length(Interval_XIC)/3));
            N_filter=N_filter_vector(P_filter(end));
            if N_filter>=10 %%%% filter length should be at most 10 scans
                N_filter=10;
            end
            Weight=ones(N_filter,1)./N_filter;
            Interval_XIC_afterfilter_1sttime=filter2(Weight,Interval_XIC); 
            Interval_XIC_afterfilter_2ndtime=filter2(Weight,Interval_XIC_afterfilter_1sttime); 
            Interval_XIC_afterfilter_3rdtime=filter2(Weight,Interval_XIC_afterfilter_2ndtime); 
            Interval_XIC_afterfilter=filter2(Weight,Interval_XIC_afterfilter_3rdtime);
%                     figure
%                     plot(intervallist{j}(1,1):intervallist{j}(1,2),Interval_XIC,'r');hold on;plot(intervallist{j}(1,1):intervallist{j}(1,2),Interval_XIC_afterfilter,'k');
%                     plot(intervallist{j}(1,1):intervallist{j}(1,2),Interval_XIC_afterfilter_2ndtime,'g');plot(intervallist{j}(1,1):intervallist{j}(1,2),Interval_XIC_afterfilter_3rdtime);
%                     grid on
%                     Height=2000;
%                     stem((intervallist_after5{j}(:,1)),2*Height*ones(length(intervallist_after5{j}(:,1)),1),'rd')
%                     stem((intervallist_after5{j}(:,2)),2*Height*ones(length(intervallist_after5{j}(:,2)),1),'kd')                   
            %%%%%
            %%%%%3. find the max peak point
            %%%%%4. find order one derivitative
            %%%%%5. find two smallest derivitative points around the max point
            
%             [max_value,max_posi, peak_start, peak_end]=Get_one_peak(Interval_XIC_afterfilter_new);
%                 intervallist_after5{j}(k,:)=[Scan_start01+peak_start-1, Scan_start01+peak_end-1]; 
%                 Maxval_posi_after5{j}(k,:)=[max_value,max_posi-peak_start];
            
            Interval_XIC_afterfilter_new=Interval_XIC_afterfilter;
            for loop_num=1:3
                [max_value(loop_num),max_posi(loop_num), peak_start(loop_num), peak_end(loop_num)]=Get_one_peak(Interval_XIC_afterfilter_new);
                Interval_XIC_afterfilter_new(peak_start(loop_num):peak_end(loop_num))=0;
            end
            ID_interval_L=find(peak_end-peak_start>=5);
            %%%%%
            if ~isempty(ID_interval_L)
                for rrr=1:length(ID_interval_L)
                    rrr1=ID_interval_L(rrr);
                    intervallist_after5{j}(num,:)=[Scan_start01+peak_start(rrr1)-1, Scan_start01+peak_end(rrr1)-1]; 
                    Maxval_posi_after5{j}(num,:)=[max_value(rrr1),max_posi(rrr1)-peak_start(rrr1)];
                    num=num+1;
                end
            end
            %%%%%

        end
    else 
        intervallist_after5{j}=[0 0];
        Maxval_posi_after5{j}=[0 0];
    end
end

IntervalList.intervallist=intervallist;%{Record(:,1)}(Record(P_m,2),:);
IntervalList.intervallist_after5=intervallist_after5;
IntervalList.Maxval_posi_after5=Maxval_posi_after5;

% colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
% figure
% for i=1:6
%     plot(xic(:,i),colorarray(i))
%     hold on
% end
% Height=1000;
% stem((intervallist_after5{1}(:,1)),2*Height*ones(length(intervallist_after5{1}(:,1)),1),'rd')
% stem((intervallist_after5{1}(:,2)),2*Height*ones(length(intervallist_after5{1}(:,2)),1),'kd')
% stem((IntervalList.intervallist_after6(:,1)),3*Height*ones(length(IntervalList.intervallist_after6(:,1)),1),'ro')
% stem((IntervalList.intervallist_after6(:,2)),3*Height*ones(length(IntervalList.intervallist_after6(:,2)),1),'ko')
% stem((IntervalList.intervallist_after7(:,1)),4*Height*ones(length(IntervalList.intervallist_after7(:,1)),1),'r+')
% stem((IntervalList.intervallist_after7(:,2)),4*Height*ones(length(IntervalList.intervallist_after7(:,2)),1),'k+')
% grid on


% IntervalList.ms2time=ms2_time;
% Exit_posi1=find(Record(:,1)==1);
% Exit_posi2=find(Record(:,1)==2);
% Exit_posi5=find(Record(:,1)==5);
% Exit_posi6=find(Record(:,1)==6);
% 
% if ~isempty(Exit_posi1) && ~isempty(Exit_posi5) && intervallist{1}(1,2)-intervallist{1}(1,1)<=size(XICs_Orbi01{1},1)/2
%     posi15=[posi15;i];
%     num15=num15+1;
% end
% 
% if ~isempty(Exit_posi2) && ~isempty(Exit_posi6) && intervallist{2}(1,2)-intervallist{2}(1,1)<=size(XICs_Orbi01{1},1)/2
%     posi26=[posi26;i];
%     num26=num26+1;
% end


clear intervallist intervallist_after5 Maxval_posi_after5            


%%%%%%%%%%%% Generate the Total Interval list based on 1st and 5th
%%%%%%%%%%%%  or 2nd and 6th in the IntervalList.intervallistv1
if iso(1)>=iso(2)

        intervallist=IntervalList.intervallist_after5;
        %%%%%%%%%%%%% raugh combine the O16 and O18 intervals
        ID_iso1=1;
        ID_iso2=5;

        Total_intervallist=[[intervallist{ID_iso1},ID_iso1*ones(size(intervallist{ID_iso1},1),1)];[intervallist{ID_iso2},ID_iso2*ones(size(intervallist{ID_iso2},1),1)]];
        Total_intervallist_new=Total_intervallist;

        IntervalList.intervallistv1=Total_intervallist_new;
else

        intervallist=IntervalList.intervallist_after5;%{Record(:,1)}(Record(P_m,2),:);
        %%%%%%%%%%%%% raugh combine the O16 and O18 intervals
        ID_iso1=2;
        ID_iso2=6;

        Total_intervallist=[[intervallist{ID_iso1},ID_iso1*ones(size(intervallist{ID_iso1},1),1)];[intervallist{ID_iso2},ID_iso2*ones(size(intervallist{ID_iso2},1),1)]];
        Total_intervallist_new=Total_intervallist;

        IntervalList.intervallistv1=Total_intervallist_new;
end

% save Totalabovev1 ID_Good_groundinterval_correction IntervalList
% load Totalabovev1


%%%%%%%%%%%%%%% interval spearman correlation correction and long's
%%%%%%%%%%%%%%% criteria
Total_interval=IntervalList.intervallistv1;
Pep_iso=iso;

for j=1:6
    XIC=XICs(:,j);
    Interval_int_threshold(j)=TOF_getThreshold_range(XIC,instrument_h_int, de_range,Times_noise_std);
end
    
ID_zero=find(Total_interval(:,1)==0);
Total_interval(ID_zero,:)=[];

if ~isempty(Total_interval)
    for k=1:size(Total_interval,1)
        Scan_start=Total_interval(k,1);
        Scan_end=Total_interval(k,2);
        ID_iso=Total_interval(k,3);
        [V_max,P_max]=max(XICs(Scan_start:Scan_end,ID_iso));
        intervaldata=zeros(Scan_end-Scan_start+1,6);
        intervalsdata01v1_afterfilter=zeros(Scan_end-Scan_start+1,6);
        N_mean_filter=5;
        for j=1:6
            intervaldata(:,j)=XICs(Scan_start:Scan_end,j);
            MS2_Peak_Intensity01(j)=intervaldata(P_max,j);

            %%%%%%%%%%%%%%%%%% XIC processed by mean filter 
            Weight=ones(N_mean_filter,1)./N_mean_filter;
            intervalsdata01v1_afterfilter(:,j)=filter2(Weight,intervaldata(:,j)); 
            MS2_Peak_Intensity01_afterfilter(j)=intervalsdata01v1_afterfilter(P_max,j);
            %%%%%%%%%%%%%%%%%%                      
        end
        
        Ms2_scannumber=P_max+Scan_start-1;
        [intervalsdata01v2 ,Scan_start_afterfilter,Scan_end_afterfilter]=MS2intervalfilter(intervaldata, MS2_Peak_Intensity01, intervalsdata01v1_afterfilter, MS2_Peak_Intensity01_afterfilter, Scan_start, Scan_end, Ms2_scannumber);
        IntervalList.intervallist_after6(k,:)=[Scan_start_afterfilter,Scan_end_afterfilter];
        %%%%%%%%%%%%%%%%%%%%%%%%%%long's peptide checks criteria

        [Judge, Labelling_efficiency]=long_criteria(Pep_iso,intervalsdata01v2,Interval_int_threshold);

        %%%%%%%%%%%%%
       IntervalList.Judge(k)=Judge;
       IntervalList.Labelling_effeciency(k)=Labelling_efficiency;
    end

    ID_long=find(IntervalList.Judge==1);
    if ~isempty(ID_long)
        IntervalList.intervallist_after7=IntervalList.intervallist_after6(ID_long,:);
        IntervalList.Labelling_effeciency_after7=IntervalList.Labelling_effeciency(ID_long);
        Jud_detect_good=1;
    else 
        IntervalList.intervallist_after7=[0 0];
        IntervalList.Labelling_effeciency_after7=0;
        Jud_detect_good=0;
    end
    
else 
    IntervalList.intervallist_after7=[0 0];
    IntervalList.Labelling_effeciency_after7=0;
    Jud_detect_good=0; 
    
end
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% combine O16 and O18 intervals based on iso mono and 1st
if Jud_detect_good==1
 
    if iso(1)>=iso(2)
        Total_intervallist=IntervalList.intervallist_after7;
        Total_intervallist_after6=IntervalList.intervallist_after6;
        %%%%%%%%%%%%% combine the O16 and O18 intervals
        ID_iso1=1;
        ID_iso2=5;
        pepXICs=XICs;
           
        [Total_intervallist_new,L]=combine_O16_O18_intervals_after7(Total_intervallist,ID_iso1,ID_iso2,pepXICs);        
        L_try=[size(Total_intervallist,1),size(Total_intervallist_new,1)];
        while L_try(end)-L_try(end-1)~=0
            [Total_intervallist_new,L]=combine_O16_O18_intervals_after7(Total_intervallist_new,ID_iso1,ID_iso2,pepXICs);
            L_try=[L_try,L];
        end
        IntervalList.intervallist_after7_combine_beforecheckL=Total_intervallist_new;
        
        [Total_intervallist_after6_new,L]=combine_O16_O18_intervals_after7(Total_intervallist_after6,ID_iso1,ID_iso2,pepXICs);        
        L_try_after6=[size(Total_intervallist_after6,1),size(Total_intervallist_after6_new,1)];
        while L_try_after6(end)-L_try_after6(end-1)~=0
            [Total_intervallist_after6_new,L]=combine_O16_O18_intervals_after7(Total_intervallist_after6_new,ID_iso1,ID_iso2,pepXICs);
            L_try_after6=[L_try_after6,L];
        end
        IntervalList.intervallist_after6_combine=Total_intervallist_after6_new;
        
    else

        Total_intervallist=IntervalList.intervallist_after7;%{Record(:,1)}(Record(P_m,2),:);
        Total_intervallist_after6=IntervalList.intervallist_after6;
        %%%%%%%%%%%%% combine the O16 and O18 intervals
        ID_iso1=2;
        ID_iso2=6;
        pepXICs=XICs;

        [Total_intervallist_new,L]=combine_O16_O18_intervals_after7(Total_intervallist,ID_iso1,ID_iso2,pepXICs);
        L_try=[size(Total_intervallist,1),size(Total_intervallist_new,1)];
        while L_try(end)-L_try(end-1)~=0
            [Total_intervallist_new,L]=combine_O16_O18_intervals_after7(Total_intervallist_new,ID_iso1,ID_iso2,pepXICs);
            L_try=[L_try,L];
        end

        IntervalList.intervallist_after7_combine_beforecheckL=Total_intervallist_new;

        [Total_intervallist_after6_new,L]=combine_O16_O18_intervals_after7(Total_intervallist_after6,ID_iso1,ID_iso2,pepXICs);        
        L_try_after6=[size(Total_intervallist_after6,1),size(Total_intervallist_after6_new,1)];
        while L_try_after6(end)-L_try_after6(end-1)~=0
            [Total_intervallist_after6_new,L]=combine_O16_O18_intervals_after7(Total_intervallist_after6_new,ID_iso1,ID_iso2,pepXICs);
            L_try_after6=[L_try_after6,L];
        end
        IntervalList.intervallist_after6_combine=Total_intervallist_after6_new;
            
    end
    clear Total_intervallist_new Total_intervallist_after6_new
else 
    IntervalList.intervallist_after7_combine_beforecheckL=[0 0];  
    IntervalList.intervallist_after6_combine=[0 0]; 
end            

%%%%%%%%%%%%%% check interval length which should be longer than 4
Id_notzeros=find(IntervalList.intervallist_after7_combine_beforecheckL(:,1)~=0);
if ~isempty(Id_notzeros)
    Id_length5=find(IntervalList.intervallist_after7_combine_beforecheckL(Id_notzeros,2)-IntervalList.intervallist_after7_combine_beforecheckL(Id_notzeros,1)>=4);
    if ~isempty(Id_length5)
        IntervalList.intervallist_after7_combine=IntervalList.intervallist_after7_combine_beforecheckL(Id_notzeros(Id_length5),:);
    else
        IntervalList.intervallist_after7_combine=[0 0];
    end
else 
    IntervalList.intervallist_after7_combine=IntervalList.intervallist_after7_combine_beforecheckL;
end
%%%%%%%%%%%%%%

%%%%%%%%%%%%%% calculate the labelling effeciency for each final intervals
minf=0.0;maxf=1;rangerate=0.6;
colorarray=['r', 'k', 'g', 'b', 'm', 'y'];
for k=1:size(IntervalList.intervallist_after7_combine,1)
    Scan_data=IntervalList.intervallist_after7_combine(k,:);
    if Scan_data(1)~=0 && Scan_data(2)~=0
        for j=1:6
            Interval_data(:,j)=xic(Scan_data(1):Scan_data(2),j);
        end
    %     figure
    %     for i=1:6
    %         plot(Interval_data(:,i),colorarray(i))
    %         hold on
    %     end
        intensity=sum(Interval_data,1);
        est=O18rateLinear(iso(1:6),intensity,minf,maxf,rangerate,1);
        IntervalList.Labelling_efficiency_after7_combine(k)=est(3);
        
    else
        IntervalList.Labelling_efficiency_after7_combine(k)=0;        
    end
    clear Interval_data
end
%%%%%%%%%%%%%%%%

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




