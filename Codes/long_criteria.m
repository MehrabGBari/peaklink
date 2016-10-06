function [Judge Labelling_efficiency]=long_criteria(iso,intervalRawData,threshold)
            Judge=0;
            isoList=iso;
            [~,maxindex]=max(isoList(1:2));
            temp0=intervalRawData(:,1);
            temp1=intervalRawData(:,2);
            temp2=intervalRawData(:,3);
            temp3=intervalRawData(:,4);
            temp4=intervalRawData(:,5);
            temp5=intervalRawData(:,6);
            temp=[temp0,temp1,temp2,temp3,temp4,temp5];
            tempmax=temp(:,maxindex);
            tempmax1=temp(:,3-maxindex);
            tempmax4=temp(:,maxindex+4);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Interval_Data=temp;            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%             for j=1:6
%                 Jud(j)=max(temp(:,j))<=threshold(j);
%             end
%             if sum(Jud)>0
%                 ID_smaller_than_th=find(Jud==1);
%                 Interval_Data(:,ID_smaller_than_th)=zeros(size(temp,1),length(ID_smaller_than_th));
%             end  
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate the labeling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% efficiency
            Spec_Pep.intervalsdata=Interval_Data;
            Spec_Pep.iso=isoList;
            %%%%%%%%%%%%%% calculate labeling eff
%             Spec_Pep_O18f=Cal_label_eff(Spec_Pep);
%             
%             I=[sum(Spec_Pep.intervalsdata,1)]';
            
            
            %%%%%%%%%%%%%%
            minf=0.0;maxf=1;rangerate=0.6;
            intensity=sum(Interval_Data,1);
            est=O18rateLinear(isoList(1:6),intensity,minf,maxf,rangerate,1);
%             est=O18rateLinear_matrix_inv(isoList(1:6),intensity,minf,maxf,rangerate,1);
%             [Spec_Pep_O18rate,Spec_Pep_O18f]=OrbitrapProduceO18RatesV1(Spec_Pep);            
            Labelling_efficiency=est(3);
            Spec_Pep_O18rate=est(2);
            
            iso_PatternO16=isoList*(1-Labelling_efficiency)^2;
            iso_PatternO18=isoList*Labelling_efficiency*(1-Labelling_efficiency)*2;
            iso_Pattern2O18=isoList*(Labelling_efficiency^2);
            iso_pattern_afterlabeling=isoList(1:6)+Spec_Pep_O18rate*(iso_PatternO16(1:6)+[0 0 iso_PatternO18(1:4)]+[0 0 0 0 iso_Pattern2O18(1:2)]);
            Normalized_Real_pattern=sum(Interval_Data,1)/sum(sum(Interval_Data,1));
            
            %%%%%%%%%%%%%%%%%%%%%% calculate the KL between the mono and
            %%%%%%%%%%%%%%%%%%%%%% 1stiso and therotical iso first 2
            %%%%%%%%%%%%%%%%%%%%%% pattern
            log_KL_Value_2iso_ther=KL_calculate(Normalized_Real_pattern(1:2),iso(1:2)/sum(iso(1:2)));
            
            if log_KL_Value_2iso_ther<=-3      
                log_KL_Value=KL_calculate(Normalized_Real_pattern,iso_pattern_afterlabeling);
                if log_KL_Value<=-3
                    Judge=1;
                else Judge=0;
                end
            else
              Judge=0;          
            end
            
            
%             Id_zeros=find(Normalized_Real_pattern~=0);
%             if ~isempty(Id_zeros)                
%                 KL_distance=sum(Normalized_Real_pattern(Id_zeros).*log(Normalized_Real_pattern(Id_zeros)./iso_pattern_afterlabeling(Id_zeros)));
%                 if log(KL_distance)<=-2.5
%                     Judge=1;
%                 else Judge=0;
%                 end
%             else 
%                 Judge=0;
%             end
    
            
            
            
            
            
            
            
%             Orbi_data01_ms2Intervalv2(i).Scanstart=Intervaldata01_Scanstart+Goodinterval_start01-1;
%             Orbi_data01_ms2Intervalv2(i).Scanend=Intervaldata01_Scanstart+Goodinterval_end01-1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%             %LC profile should be simliar for first and fifth XIC interval,
%             %first and second; 5th and 6th
%             corr1=cosinecorr(tempmax,tempmax1);
%             corr2=cosinecorr(temp(:,maxindex+4),temp(:,7-maxindex));
%             corr3=cosinecorr(temp(:,maxindex),temp(:,maxindex+4));
%             tcount=0;
%             for nn=1:size(temp,1)
%                 if KLdistance(temp0(nn),temp1(nn),isoList(1),isoList(2))<-2.5 && KLdistance(temp2(nn), temp3(nn),sum(temp2), sum(temp3))<-2.5 && KLdistance(temp4(nn),temp5(nn),sum(temp4), sum(temp5))<-2.5 &&  KLdistance(tempmax(nn),tempmax4(nn),sum(tempmax), sum(tempmax4))<-2.5
%                     tcount=tcount+1;
%                 end
%             end
% 
%             if corr1^2>threshold && corr2^2>threshold && corr3^2>threshold && KLdistance(isoList(1),isoList(2),sum(temp0), sum(temp1))<-2.5
%                 if tcount>4 || 2*tcount>size(temp,1)
%                     Judge=1;
%                 end
%             end
%             %%%%%%%%%%%%%