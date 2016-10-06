function  [max_value,max_posi, peak_start, peak_end]=Get_one_peak(Interval_XIC_afterfilter)
                    
%                     figure
%                     plot(Interval_XIC_afterfilter)
                    [V_max_aftfilter,P_max_aftfilter]=max(Interval_XIC_afterfilter);
                    Derivat_value=diff(Interval_XIC_afterfilter);
                    if P_max_aftfilter~=1                        
                        Po_forward_value=find(Derivat_value<=0);
                        ID_forward_min=find(Po_forward_value<=(P_max_aftfilter-1));
                        if isempty(ID_forward_min)
                            peak_start=1;
                        else
                            peak_start=Po_forward_value(ID_forward_min(end))+1;
                        end
                    else peak_start=1;
                    end
                    if P_max_aftfilter~=length(Interval_XIC_afterfilter)                        
                        Po_backward_value=find(Derivat_value>=0);
                        ID_backward_min=find(Po_backward_value>(P_max_aftfilter-1));
                        if isempty(ID_backward_min)
                            peak_end=length(Interval_XIC_afterfilter);
                        else
                            peak_end=Po_backward_value(ID_backward_min(1));
                        end
                    else peak_end=length(Interval_XIC_afterfilter);
                    end
                    
                    max_value=V_max_aftfilter;
                    max_posi=P_max_aftfilter;
                    
                    
                    