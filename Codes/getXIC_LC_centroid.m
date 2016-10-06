function XICs=getXIC_LC_centroid(peakl,centermz,tolerance)
% we expect windows is a length N*2 matrix representing the start and end
% of N windows
%t=cputime;

windowNum=length(centermz);
%overlapNum=5;
XICs=zeros(length(peakl.scan),windowNum);

massdiff=tolerance*centermz*1e-6;
windowleft=centermz-massdiff;
windowright=centermz+massdiff;

numoverlap=zeros(size(centermz));
for i=1:length(centermz)-1
    for j=i+1:length(centermz)
        if windowleft(j)<=windowright(i) 
            numoverlap(i)=numoverlap(i)+1;
        else
            break;
        end
    end
end
%matlabpool;
for scanNumber = 1 : length(peakl.scan)
    mz = peakl.scan(scanNumber).peak(:,1);
    int = peakl.scan(scanNumber).peak(:,2);
    windowidx=1;
    tempxic=zeros(1,windowNum);
    for mzidx=1:length(mz)
    	currentmz=mz(mzidx);
    	while(currentmz>=windowright(windowidx))
    		windowidx=windowidx+1;
            if windowidx>windowNum
                break;
            end
    	end
    	if windowidx>windowNum
    		break;
    	end
    	if(currentmz>windowleft(windowidx))
    		%XICs(scanNumber,windowidx)=XICs(scanNumber,windowidx)+int(mzidx);
            tempxic(windowidx)=max(tempxic(windowidx),int(mzidx));
            if numoverlap(windowidx)>0
                for i=1:numoverlap(windowidx)
                    if (windowidx+i<=windowNum)
                        if(currentmz>windowleft(windowidx+i))
                            %XICs(scanNumber,windowidx+i)=XICs(scanNumber,windowidx+i)+int(mzidx);
                            tempxic(windowidx+i)=max(tempxic(windowidx+i),int(mzidx));
                        end
                    end
                end
            end
    	end
    end
    XICs(scanNumber,:)=tempxic;
end    
%matlabpool close;
%cputime-t
    