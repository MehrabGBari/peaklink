function [TestAccuracy_SVM]=BuildAWModelParameters_1_test(Testing_Information_Matrix,Model_Testing,parameters,SVMStruct)


% % % %  %%%%%%%%%%%%%%%%%% file count 1,2 and 3 is for model AB AC BC
% % % %  %%%%%%%%%%%%%%%%%%
% % % % 
% % % %       Corresponding_AW=[];Corresponding_PATc=[];Corresponding_PATnc=[];
% % % %       NonCorresponding_AW=[];NonCorresponding_PATnc=[];NonCorresponding_PATc=[];
% % % %     
% % % %     for i=1:length(Model_Training)
% % % %       
% % % %         MS2_ID01=Training_Information_Matrix{i}{1}.MS2_ID;
% % % %         MS2_ID02=Training_Information_Matrix{i}{2}.MS2_ID;   
% % % % 
% % % %         AW_Matrix=Model_Training{i}.AW;
% % % %         PATc_Matrix=Model_Training{i}.PAT_c; 
% % % %         PATnc_Matrix=Model_Training{i}.PAT_nc;
% % % %         
% % % %         Corresponding_AW(i)=AW_Matrix(1,MS2_ID02);
% % % %         Corresponding_PATc(i)=PATc_Matrix(1,MS2_ID02);%PAT corr using corr param
% % % %         Corresponding_PATnc(i)=PATnc_Matrix(1,MS2_ID02);%PAT corr using Noncorr param
% % % %         %V1=AT_Matrix(:,MS2_ID02);V1(MS2_ID01)=[];
% % % %         V2=AW_Matrix;V2(MS2_ID02)=[];
% % % %         NON_Zero=V2~=0;
% % % %         V2=V2(NON_Zero);
% % % %         NonCorresponding_AW=[NonCorresponding_AW, V2];
% % % %         
% % % %         V3=PATnc_Matrix;V3(MS2_ID02)=[];
% % % %         NON_Zero1=V3~=0;
% % % %         V3=V3(NON_Zero1);
% % % %         NonCorresponding_PATnc=[NonCorresponding_PATnc, V3];%PAT noncorr using Noncorr param
% % % %         
% % % %         V4=PATc_Matrix;V4(MS2_ID02)=[];
% % % %         NON_Zero2=V4~=0;
% % % %         V4=V4(NON_Zero2);
% % % %         NonCorresponding_PATc=[NonCorresponding_PATc, V4];%PAT noncorr using corr param
% % % %                
% % % %     end
% % % %     NON_ZeroAW=Corresponding_AW~=0;%removing zero amunts
% % % %     Corresponding_AW=Corresponding_AW(NON_ZeroAW);
% % % %     Corresponding_PATc=Corresponding_PATc(NON_ZeroAW);
% % % %     Corresponding_PATnc=Corresponding_PATnc(NON_ZeroAW);
% % % %     
% % % %     AW_PAWMHAT_sig=gamfit(Corresponding_AW);
% % % %     AW_PAWMHAT_notsig=gamfit(NonCorresponding_AW);

    %%%parameters.PAWa_AWmodel=[AW_PAWMHAT_sig,AW_PAWMHAT_notsig];

   
    for i=1:length(Model_Testing)
      
        for j=1:length(Model_Testing{i}.AW)
            if Model_Testing{i}.AW(j)~=0
                PAW_corr_param (j)= pdf('Gamma',Model_Testing{i}.AW(j),parameters.PAWa_AWmodel(2),parameters.PAWa_AWmodel(1));
                PAW_noncorr_param (j)= pdf('Gamma',Model_Testing{i}.AW(j),parameters.PAWa_AWmodel(4),parameters.PAWa_AWmodel(3));
                PAWc_div_PAWnc(j)= PAW_corr_param (j)/ PAW_noncorr_param (j);
                PAWnc_div_PAWc(j)= PAW_noncorr_param (j)/ PAW_corr_param (j);
                PATc_div_PATnc(j)=Model_Testing{i}.PAT_c(j)/Model_Testing{i}.PAT_nc(j);
            else
                PAW_corr_param (j)= 0;
                PAW_noncorr_param (j)= 0;
                PAWc_div_PAWnc(j)= 0;
                PAWnc_div_PAWc(j)= 0;
                PATc_div_PATnc(j)=0;
            end
        end

                 Result{i}.AWc=PAW_corr_param/max(PAW_corr_param);
                 Result{i}.AWnc=PAW_noncorr_param/max(PAW_noncorr_param);
                 if max(PAWc_div_PAWnc)~=0
                 Result{i}.PAWc_div_PAWnc=PAWc_div_PAWnc/max(PAWc_div_PAWnc);
                 else
                     Result{i}.PAWc_div_PAWnc=PAWc_div_PAWnc;
                 end
                 Result{i}.PAWnc_div_PAWc=PAWnc_div_PAWc/max(PAWnc_div_PAWc);
                 Result{i}.PATc_div_PATnc=PATc_div_PATnc;%/max(PATc_div_PATnc);
          clear PAW_corr_param PAW_noncorr_param PAWc_div_PAWnc PAWnc_div_PAWc PATc_div_PATnc
    end
    
      Corresponding_PAW=[];Corresponding_PAT=[];
      NonCorresponding_PAW=[];NonCorresponding_PAT=[];
    
        for i=1:length(Result)
      
        MS2_ID01=Testing_Information_Matrix{i}{1}.MS2_ID;
        MS2_ID02=Testing_Information_Matrix{i}{2}.MS2_ID;   

        PAW_Matrix=Result{i}.PAWc_div_PAWnc;
        PAT_Matrix=Result{i}.PATc_div_PATnc;
        
        Corresponding_PAW(i)=PAW_Matrix(1,MS2_ID02);
        Corresponding_PAT(i)=PAT_Matrix(1,MS2_ID02);%PAT corr using corr param
        
        %V1=AT_Matrix(:,MS2_ID02);V1(MS2_ID01)=[];
        V2=PAW_Matrix;V2(MS2_ID02)=[];
        NON_Zero=V2~=0;
        V2=V2(NON_Zero); 
        NonCorresponding_PAW=[NonCorresponding_PAW, V2];
        
        V3=PAT_Matrix;V3(MS2_ID02)=[];
        NON_Zero1=V3~=0;
        V3=V3(NON_Zero1);
        NonCorresponding_PAT=[NonCorresponding_PAT, V3];%PAT noncorr using Noncorr param
        
         
        end
%          NON_ZeroPAW=Corresponding_PAW~=0;%removing zero amunts
%     Corresponding_PAW=Corresponding_PAW(NON_ZeroPAW);
%     Corresponding_PAT=Corresponding_PAT(NON_ZeroPAW);
%  NonCorresponding_PAT=(NonCorresponding_PAT-mean(NonCorresponding_PAT))/std(NonCorresponding_PAT);          
% Corresponding_PAT=(Corresponding_PAT-mean(Corresponding_PAT))/std(Corresponding_PAT);  
% NonCorresponding_PAW=(NonCorresponding_PAW-mean(NonCorresponding_PAW))/std(NonCorresponding_PAW);  
% Corresponding_PAW=(Corresponding_PAW-mean(Corresponding_PAW))/std(Corresponding_PAW);  

  L1=min(length(NonCorresponding_PAT),length(NonCorresponding_PAW));
  L2=min(length(Corresponding_PAT),length(Corresponding_PAW));
%figure;plot(NonCorresponding_PAT(1:L1),NonCorresponding_PAW(1:L1),'bo',Corresponding_PAT(1:L2),Corresponding_PAW(1:L2),'r*');title('V: PAWc/PAWnc , H: PATc/PATnc');axis([-.1 9 -0.1 1.1])
TestData=[Corresponding_PAT(1:L2) NonCorresponding_PAT(1:L1);Corresponding_PAW(1:L2) NonCorresponding_PAW(1:L1)];
N1=length(TestData);n_train=length(Corresponding_PAW);
T=[ones(n_train,1);zeros(N1-n_train,1)]; % the label vector

%figure;SVMStruct = svmtrain(TrainData,T,'KERNEL_FUNCTION','rbf','BOXCONSTRAINT',30,'showplot',true);
GROUP =  svmclassify(SVMStruct,TestData');
%TestAccuracy_SVM=mean(double(GROUP==T)*100)

TestAccuracy_SVM=mean(double(GROUP(1:L2,1)==T(1:L2,1))*100)

