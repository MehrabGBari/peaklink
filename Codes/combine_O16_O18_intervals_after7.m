function   [Total_intervallist_new,L]=combine_O16_O18_intervals_after7(Total_intervallist,ID_iso1,ID_iso2,pepXICs)
            
            [V,P]=sort(Total_intervallist(:,1));
            Total_intervallist_sort=Total_intervallist(P,:);
            Posi_record=ones(1,size(Total_intervallist_sort,1));
            for i=1:size(Total_intervallist_sort,1)-1

                Indicator=find(Total_intervallist_sort(i+1:end,1)<=Total_intervallist_sort(i,2));

                if ~isempty(Indicator)
                    Indicator=Indicator+i;
                    for id=1:length(Indicator)
                        L1=Total_intervallist_sort(i,2)-Total_intervallist_sort(i,1)+1;
                        L2=Total_intervallist_sort(Indicator(id),2)-Total_intervallist_sort(Indicator(id),1)+1;
                        L3=min([Total_intervallist_sort(i,2),Total_intervallist_sort(Indicator(id),2)])-Total_intervallist_sort(Indicator(id),1)+1;
%                         if Total_intervallist_sort(Indicator(id),2)<=Total_intervallist_sort(i,2)
                        if L3/L1>=0.8 || L3/L2>=0.8
                            if L1>=L2
                                Posi_record(Indicator(id))=0;
                            else
                                Posi_record(i)=0;
                            end
                        end
                    end
                end
  
            end
            Total_intervallist_sort(Posi_record==0,:)=[];
            
            
%             Total_Maxval_posi_after5_sort=Total_Maxval_posi_after5(P,:);
%             Max_Peak_posi=Total_intervallist_sort(:,1)+Total_Maxval_posi_after5_sort(:,2);
            com_Peak_length=Total_intervallist_sort(2:end,2)-Total_intervallist_sort(1:end-1,1);
            Interval_overlap=Total_intervallist_sort(2:end,1)-Total_intervallist_sort(1:end-1,2);
            ID_Interval_overlap=find(Interval_overlap<0);
%             Max_peak_diff=Max_Peak_posi(2:end,1)-Max_Peak_posi(1:end-1,1);
%             ID_Interval_overlap=find(Max_peak_diff<5);

            if ~isempty(ID_Interval_overlap)
                ID_delete=[];
                for k=1:length(ID_Interval_overlap) 
                    Thres=min([com_Peak_length(ID_Interval_overlap(k))*0.3,5]);
                    
                    Scan_forward=Total_intervallist_sort(ID_Interval_overlap(k),:);
                    Scan_backward=Total_intervallist_sort(ID_Interval_overlap(k)+1,:);

                    Interval_XIC_mono=pepXICs(Scan_forward(1):Scan_forward(2),ID_iso1);
                    Interval_XIC_5th=pepXICs(Scan_backward(1):Scan_backward(2),ID_iso2);
                    
                    N_filter=5;
                    Weight=ones(N_filter,1)./N_filter;
                    Interval_XIC_mono_afterfilter=filter2(Weight,Interval_XIC_mono); 
                    Interval_XIC_5th_afterfilter=filter2(Weight,Interval_XIC_5th); 
                    
                    [v_max_iso1,P_max_iso1]=max(Interval_XIC_mono_afterfilter);
                    [v_max_iso2,P_max_iso2]=max(Interval_XIC_5th_afterfilter);
                    
                    Peak_diff=abs(Scan_forward(1)+P_max_iso1-1-Scan_backward(1)-P_max_iso2+1);
                    if Peak_diff<Thres                    
                        if max(Interval_XIC_mono)>=max(Interval_XIC_5th)
                                ID_delete=[ID_delete;ID_Interval_overlap(k)+1];
                        else ID_delete=[ID_delete;ID_Interval_overlap(k)];
                        end
                    end
                end
                Total_intervallist_new=Total_intervallist_sort;
                Total_intervallist_new(ID_delete,:)=[];
            else  Total_intervallist_new=Total_intervallist_sort;
            end
            
            
            for i=1:size(Total_intervallist_new,1)-1

                Indicator=(Total_intervallist_new(i+1,1)<=Total_intervallist_new(i,2));

                if  Indicator
     
                    New_scan=floor((Total_intervallist_new(i,2)+Total_intervallist_new(i+1,1))/2);
                    Total_intervallist_new(i,2)=New_scan;
                    Total_intervallist_new(i+1,1)=New_scan+1;



                end
  
            end

            
%             Total_intervallist_new;
            L=size(Total_intervallist_new,1);
            
%             Interval_overlap_2nd=Total_intervallist_new(2:end,1)-Total_intervallist_new(1:end-1,2);
%             ID_over=find(Interval_overlap_2nd<0);
%             if ~isempty(ID_over)
%                 ID_delete_2nd=[];
%                 for k=1:length(ID_over)                    
%                         Scan_forward_2nd=Total_intervallist_new(ID_Interval_overlap(k),:);
%                         Scan_backward_2nd=Total_intervallist_new(ID_Interval_overlap(k)+1,:);
%                     
%                         Interval_XIC_mono_2nd=pepXICs(Scan_forward_2nd(1):Scan_forward_2nd(2),1);
%                         Interval_XIC_5th_2nd=pepXICs(Scan_backward_2nd(1):Scan_backward_2nd(2),5);
%                         if max(Interval_XIC_mono_2nd)>=max(Interval_XIC_5th_2nd)
%                                 ID_delete_2nd=[ID_delete_2nd;ID_over(k)+1];
%                         else ID_delete_2nd=[ID_delete_2nd;ID_over(k)];
%                         end
%                 end
%                 Total_intervallist_new_2nd=Total_intervallist_new;
%                 Total_intervallist_new_2nd(ID_delete,:)=[];
%             else  Total_intervallist_new_2nd=Total_intervallist_new;
%             end
            

