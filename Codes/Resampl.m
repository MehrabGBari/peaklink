function ReSample_Sig=Resampl(Signal)
if nnz(Signal)<2
    ReSample_Sig=Signal;%do not perform resampling for signal length<2
else
   [nr,~]=size(Signal);
     TimeA = 0:nr-1;
     TimeA1 = 0:.25:nr; 
ReSample_Sig=zeros(length(TimeA1)-4,1);
%for IntsIndex=1
     Signal1=Signal;%(:,IntsIndex);
%      if ~(sum(Signal1)==0)
%      Signal1=Signal1/sum(Signal1); 
%      end
     ReSample_Sig1= interp1(TimeA,Signal1,TimeA1);ReSample_Sig1=ReSample_Sig1(1,1:end-4)';
ReSample_Sig=ReSample_Sig1; 
%end

end            