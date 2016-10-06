function [newpeplist, newchargestate, groundtruthinterval, IntervalList, ms2time, newmonoXICs,newiso1stXICs]=verifypeptide(pep,chargestate,ms2information,monoXICs,iso1stXICs,retentiontl1)

mininterval=4;
maxnointervals=6;

% ms2time=ms2information;
posi=[];
num=1;
for i=1:size(monoXICs,2)
        ms2_time=ms2information(i);
        xic=monoXICs(:,i);
        [intervallist]=intervaldetection(xic,mininterval,maxnointervals);

        if intervallist(1,1)~=0 && intervallist(1,2)~=0
                for j=1:size(intervallist,1)
                    start02=retentiontl1(intervallist(j,1));
                    end02=retentiontl1(intervallist(j,2));
                    if ms2_time>=start02 && ms2_time<=end02
                        posi=[posi;i, j];    %%%% i for peptide;  j for interval
                        ms2time(num)=ms2_time;
                        IntervalList{num}.intervallist=intervallist;
                        num=num+1;
                    end
                end
        end
        
        
end

newpeplist=pep(posi(:,1));
newmonoXICs=monoXICs(:,posi(:,1));
newiso1stXICs=iso1stXICs(:,posi(:,1));
% newiso2ndXICs=iso2ndXICs(:,posi(:,1));
newchargestate=chargestate(posi(:,1));
groundtruthinterval=posi;




































