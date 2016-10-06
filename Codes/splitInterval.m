function [intervalList,intervalCount]=splitInterval(smoothXIC,diffXIC,intervalList,intervalCount, minLCLength)
% Since each interval is split to at most two intervals create a temp
% intervallist 
newIntervalList=zeros(intervalCount*10,2);
newIntervalIndex=0;
for intervalIndex=1:intervalCount
    intervalLength=intervalList(intervalIndex,2)-intervalList(intervalIndex,1)+1;
    zeroCrossingPoints=ones(1,intervalLength+2);
    zCPheight=ones(1,intervalLength+2);
    zCPid=1;
    zeroCrossingPoints(1)=intervalList(intervalIndex,1);
    zCPheight(1)=smoothXIC(intervalList(intervalIndex,1));
    for tempid=intervalList(intervalIndex,1)+1:intervalList(intervalIndex,2)-1
        if (sign(diffXIC(tempid))+sign(diffXIC(tempid-1)))==0 || diffXIC(tempid)==0
            
            zeroCrossingPoints(zCPid)=tempid;
            zCPheight(zCPid)=smoothXIC(tempid);
            zCPid=zCPid+1;
        end   
    end
    zeroCrossingPoints(zCPid)=intervalList(intervalIndex,2);
    zCPheight(zCPid)=smoothXIC(intervalList(intervalIndex,2));
    zeroCrossingPoints=zeroCrossingPoints(1:zCPid);
    zCPheight=zCPheight(1:zCPid);
    valleyPoints=zeros(1,zCPid);
    valleyPointCount=0;
    if zCPid>1
        % find all valley points
        for count=2:zCPid-1
            if zCPheight(count)< zCPheight(count-1) && zCPheight(count)< zCPheight(count+1)
                valleyPointCount=valleyPointCount+1;
                valleyPoints(valleyPointCount)=zeroCrossingPoints(count);
            end
        end
        if valleyPointCount>0
            % first replace the origial interval end point with the first valley
            % point
            oldIntervalEnd=intervalList(intervalIndex,2);
           intervalList(intervalIndex,2)=valleyPoints(1);
           % next generate a new interval
           newIntervalIndex=newIntervalIndex+1;
           % Set the start point of the new interval
           newIntervalList(newIntervalIndex,1)=valleyPoints(1)+1;
           if valleyPointCount>1 
                % generate a new intervallist based on the valley points. 
                for vpindex=2:valleyPointCount
                    newIntervalList(newIntervalIndex,2)=valleyPoints(vpindex);
                    newIntervalIndex=newIntervalIndex+1;
                    newIntervalList(newIntervalIndex,1)=valleyPoints(vpindex)+1;
                end
           end 
           newIntervalList(newIntervalIndex,2)=oldIntervalEnd;
        end
    end
end
newIntervalList=newIntervalList(1:newIntervalIndex,:);
intervalList=[intervalList;  newIntervalList];
intervalCount=intervalCount+newIntervalIndex;

% Get rid of intervals with short length


newintervalList=zeros(intervalCount,2);
newintervalCount=0;
 for index2=1:intervalCount
     intervalLength=intervalList(index2,2)-intervalList(index2,1)+1;
     if intervalLength>minLCLength;
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
 
 heightv=zeros(intervalCount,1);
 for i=1:intervalCount
     heightv(i)=max(smoothXIC(intervalList(i,1):intervalList(i,2)));
 end    
 
 %intervallength=intervallength(iids);
 %[dummy, sortid]=sort(intervallength,'descend');
 [dummy, sortid]=sort(heightv,'descend');
 intervalList=intervalList(sortid,:);
  
        
        
       
    
    
        
    
