function [Judge Labelling_efficiency]=checkDataModelO18(isoList,intervalRawData, klthreshold)        
           Labelling_efficiency=0;
           isoList=isoList(1:6)./sum(isoList(1:6));
           
           % Estimate the abundance at the three positions using modified
           % Yao's method.
           
           [AbundanceO16O18 Labelling_efficiency]=estimateO16O18ModifiedYao(isoList,intervalRawData);
           if Labelling_efficiency>0
               AHeavy=AbundanceO16O18(3)/Labelling_efficiency^2;
               if (AbundanceO16O18(1)<AHeavy*(1-Labelling_efficiency)^2)
                   AbundanceO16O18(1)=AHeavy*(1-Labelling_efficiency)^2;
               end  
           end
            collapsedScan=sum(intervalRawData,1);
            meanOb=mean(collapsedScan);
            stdOb=std(collapsedScan);
            if 2*stdOb/meanOb <abs(isoList(1)-isoList(2)) % the peaks have equal height, then treat it as noise
                Judge=0;
                Labelling_efficiency=0;
                return;
            end    
            
            if Labelling_efficiency<0.67
                Judge=0;
                Labelling_efficiency=0;
                return; 
            end
           
           estimatedPattern1=AbundanceO16O18(1)*isoList(1:6);
           estimatedPattern2=[0 0 AbundanceO16O18(2)*isoList(1:4)]; 
           estimatedPattern3=[0 0 0 0 AbundanceO16O18(3)*isoList(1:2)]; 
           estimatedPattern=estimatedPattern1+estimatedPattern2+estimatedPattern3;
           
          
           estimatedPattern=estimatedPattern./sum(estimatedPattern);
           obervedPattern=collapsedScan./sum(collapsedScan);
           klscore= KL_calculate(estimatedPattern,obervedPattern);
           % simulate random patterns with equal variance to that of the
           % observed pattern and mean, check the mean kl score 
           
           if klscore<klthreshold 
              Judge=1;
             
           else 
              Judge=0;
           end   
           
         
            
%                 %%%%%%%%%%%%%%
%                 minf=0.0;maxf=1;rangerate=0.6;
%                 intensity=sum(Interval_Data,1);
%                 est=O18rateLinear(isoList(1:6),intensity,minf,maxf,rangerate,1);
%     %             est=O18rateLinear_matrix_inv(isoList(1:6),intensity,minf,maxf,rangerate,1);
%     %             [Spec_Pep_O18rate,Spec_Pep_O18f]=OrbitrapProduceO18RatesV1(Spec_Pep);            
%                 Labelling_efficiency=est(3);
%                 Spec_Pep_O18rate=est(2);
% 
%                 iso_PatternO16=isoList*(1-Labelling_efficiency)^2;
%                 iso_PatternO18=isoList*Labelling_efficiency*(1-Labelling_efficiency)*2;
%                 iso_Pattern2O18=isoList*(Labelling_efficiency^2);
%                 iso_pattern_afterlabeling=isoList(1:6)+Spec_Pep_O18rate*(iso_PatternO16(1:6)+[0 0 iso_PatternO18(1:4)]+[0 0 0 0 iso_Pattern2O18(1:2)]);
%                 Normalized_Real_pattern=sum(Interval_Data,1)/sum(sum(Interval_Data,1));
% 
%                 %%%%%%%%%%%%%%%%%%%%%% calculate the KL between the mono and
%                 %%%%%%%%%%%%%%%%%%%%%% 1stiso and therotical iso first 2
%                 %%%%%%%%%%%%%%%%%%%%%% pattern
%                 log_KL_Value_2iso_ther=KL_calculate(Normalized_Real_pattern(1:2),iso(1:2)/sum(iso(1:2)));
% 
%                 if log_KL_Value_2iso_ther<=-3      
%                     log_KL_Value=KL_calculate(Normalized_Real_pattern,iso_pattern_afterlabeling);
%                     if log_KL_Value<=-3
%                         Judge=1;
%                     else Judge=0;
%                     end
%                 else
%                   Judge=0;          
%                 end
% 
% 
%     %             Id_zeros=find(Normalized_Real_pattern~=0);
%     %             if ~isempty(Id_zeros)                
%     %                 KL_distance=sum(Normalized_Real_pattern(Id_zeros).*log(Normalized_Real_pattern(Id_zeros)./iso_pattern_afterlabeling(Id_zeros)));
%     %                 if log(KL_distance)<=-2.5
%     %                     Judge=1;
%     %                 else Judge=0;
%     %                 end
%     %             else 
%     %                 Judge=0;
%     %             end
% 

            
            
            
            
            
            
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