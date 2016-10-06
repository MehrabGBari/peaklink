function [Corr_all Corr_9 Corr_7]=wavedec_corr(SignalA,SignalB)
% wavelet decomposition and normalizing signal by sum sig+coefftiont correlation 
%  SignalA=alignedProfile_A_Corr;
%  SignalB=alignedProfile_B_Corr;
Ra=nnz(SignalA);
Rb=nnz(SignalB);

         switch Ra*Rb
                    case 0 
                        if (Ra==0 && Rb==0)
                            Corr_all=1;Corr_9=1;Corr_7=1;
                        elseif (Ra==1 || Rb==1)
                            Corr_all=0.5;Corr_9=0.5;Corr_7=0.5;
                        else
                         Corr_all=0;Corr_9=0;Corr_7=0;
                        end
                    case 1
                         Corr_all=1;Corr_9=1;Corr_7=1;
             otherwise
               if ((Ra==1 || Rb==1) && (abs(Ra-Rb)>=5))
                    Corr_all=0;Corr_9=0;Corr_7=0;
               else
                    SignalA=SignalA/sum(SignalA);  
                    SignalB=SignalB/sum(SignalB); 
                    
                    [CA1,La] = wavedec(SignalA,6,'db16');
                    [CB1,Lb] = wavedec(SignalB,6,'db16');
                               END=min(length(CA1),length(CB1));
                    
                     Corr_all=abs(corr(CA1(1:END,1),CB1(1:END,1)));
%                      Corr_13=abs(corr(CA1(1:12,1),CB1(1:12,1)));
%                      Corr_11=abs(corr(CA1(1:10,1),CB1(1:10,1)));
                     Corr_9=abs(corr(CA1(1:9,1),CB1(1:9,1)));
                     Corr_7=abs(corr(CA1(1:7,1),CB1(1:7,1)));
%                      Corr_5=abs(corr(CA1(1:8,1),CB1(1:8,1)));
%                      Corr_3=abs(corr(CA1(1:6,1),CB1(1:6,1)));
%                      Corr_1=abs(corr(CA1(1:1,1),CB1(1:1,1)));
%                     
               end
        end
 
  
