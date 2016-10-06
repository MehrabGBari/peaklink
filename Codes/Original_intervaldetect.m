function IntervalList=Original_intervaldetect(Orbi_peplist01_final,Orbi_ms2information01_final,Orbi_monoXICs01,Orbi_retentiont01l1)
mininterval=4;
maxnointervals=20;
posi=[];
num=0;
for i=1:length(Orbi_peplist01_final)
        ms2_time=Orbi_ms2information01_final(i);        
        Record=[];
        MaxValue=[];
        xic=zeros(size(Orbi_monoXICs01,1),6);
        xic_afterfilter=xic;
        for j=1%:6
            xic(:,j)=Orbi_monoXICs01(:,i);
%             %%%%%%%%%%%%%%%%%% XIC processed by mean filter 
%             Weight=[0.2 0.2 0.2 0.2 0.2]';
%             xic_afterfilter(:,j)=filter2(Weight,xic(:,j)); 
%             intervallist{j}=intervaldetection(xic_afterfilter(:,j),mininterval,maxnointervals);
%             %%%%%%%%%%%%%%%%%%
            intervallist{j}=intervaldetection(xic(:,j),mininterval,maxnointervals);
            if intervallist{j}(1,1)~=0 && intervallist{j}(1,2)~=0
                for k=1:size(intervallist{j},1)
                    Scan_start01=intervallist{j}(k,1);
                    Scan_end01=intervallist{j}(k,2);
                    if Scan_end01-Scan_start01<=length(Orbi_monoXICs01(:,i))*0.75
                        start01=Orbi_retentiont01l1(Scan_start01);
                        end01=Orbi_retentiont01l1(Scan_end01);
                        if ms2_time>=start01 && ms2_time<=end01
                            IntervalList(i,:)=[Scan_start01,Scan_end01];
                        end
                    end
                end
            end
        end
               
%         if ~isempty(MaxValue)
%             [V_m,P_m]=max(MaxValue);
%             IntervalList{i}.intervallist=intervallist;%{Record(:,1)}(Record(P_m,2),:);
%             IntervalList{i}.Record=Record;
%             IntervalList{i}.MaxValue=MaxValue;
%             IntervalList{i}.ms2time=ms2_time;
%             Exit_posi1=find(Record(:,1)==1);
%             Exit_posi4=find(Record(:,1)==4);
%             Exit_posi2=find(Record(:,1)==2);
%             Exit_posi5=find(Record(:,1)==5);
%             if (~isempty(Exit_posi1) && ~isempty(Exit_posi4)) || (~isempty(Exit_posi2) && ~isempty(Exit_posi5))
%                 posi=[posi;i];
%                 num=num+1;
%             end
%         end
            
end

