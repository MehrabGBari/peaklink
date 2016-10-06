% Returns a list a intervals of the xic ordered by peak height
%Inputs: XIC, a column sigle array;
%th: threshold to be applied to xic;
%indicate maximum gap that can be ignored between intervals of ten.
%mininterval: Mimimu length of an interval.




function [intervalList, intervalCount]=getInterval(xic,th,mininterval)
ll=length(xic);
intervalmap=zeros(ll,1);
if isempty(th)==1
    th=20;
end    

for i=1:ll
    if xic(i)>th;
        intervalmap(i)=1;
    end
end    
%ids=find(xic>=th);
%intervalList=zeros(maxnointervals,2);
if sum(intervalmap)==0
    intervalList=[0 0];
    intervalCount=0;
    return;
else    
intervalList=zeros(ll,2);
 intervalCount=0;

 % get the difference in the map intervalmap(i)-intervalmap(i-1)
 diff=getDiff(intervalmap,ll,-1);

 for index=1:ll
     if diff(index)==1
         intervalCount=intervalCount+1;
         intervalList(intervalCount,1)=index;
     end
     if diff(index)==-1
         intervalList(intervalCount,2)=index-1;
     end    
     
 end    
 if intervalCount>0
      intervalList=intervalList(1:intervalCount,:);
 else 
     intervalList=[0 0];
     intervalCount=0;
     return;
 end
 newintervalList=zeros(intervalCount,2);
 newintervalCount=0;
 for index2=1:intervalCount
     intervalLength=intervalList(index2,2)-intervalList(index2,1)+1;
     if intervalLength>mininterval;
         newintervalCount=newintervalCount+1;
         newintervalList(newintervalCount,:)=intervalList(index2,:);
     end
 end
 if newintervalCount>0
  intervalList=newintervalList(1:newintervalCount,:);
  intervalCount=newintervalCount;
 else
     intervalList=[0 0];
     intervalCount=0;
     return;
 end    
 
% heightv=zeros(intervalCount,1);
% for i=1:intervalCount
%     heightv(i)=max(xic(intervalList(i,1):intervalList(i,2)));
% end    
 
 %intervallength=intervallength(iids);
 %[dummy, sortid]=sort(intervallength,'descend');
 %[dummy, sortid]=sort(heightv,'descend');
 %intervalList=intervalList(sortid,:);
  
 end   
end
 
 
 
 
