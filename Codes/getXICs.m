% Version 1.1 Created by Jianqiu Zhang Feb. 2010
% Extract many XICs given an array of centermz values and tolerance in ppm.
% Input parameters: 
% datascan- the datastructure of LC/MS scans.
% centermz- an array of m/z values
% tolerance- the tolerance for mass difference of the instrument in ppm.
% Output parameters: 
% XIC an matrix with size S*N, where S is the total number of scans, and N
% is the length of centermz.
  
function XICs=getXICs(datascan,centermz,tolerance)


N=length(centermz);
XICs=zeros(length(datascan.scan),N);
massdiff=tolerance*centermz*1e-6;
for scanNumber = 1 : length(datascan.scan)
%     mz = datascan.scan(scanNumber).peaks.mz(1 : 2 : end - 1);
%     int = datascan.scan(scanNumber).peaks.mz(2 : 2 : end);
    mz = datascan.scan(scanNumber).peak(:,1);
    int = datascan.scan(scanNumber).peak(:,2);
    for w=1:N
        ids=find(mz>(centermz(w)-massdiff(w)) & mz<(centermz(w)+massdiff(w)));
        if isempty(ids)==0
           %[minmz minmzid]=min(abs(mz(ids)-centermz(windowid)));
           %inttemp=int(ids);
           XICs(scanNumber,w)=sum(int(ids));
        else   
          XICs(scanNumber,w)=0;
       end   
    end   
end    
    