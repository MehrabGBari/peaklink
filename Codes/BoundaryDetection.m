function [intervalsdata01v2 ,Scan_start_afterfilter,Scan_end_afterfilter]=BoundaryDetection(intervaldata, Scan_start,Ms2_scannumber,klthreshold)
    
   % first get template at the top peak position given by the scannumber
   maxposi=Ms2_scannumber-Scan_start+1;
   template=intervaldata(maxposi,:);
   IntervalSize=size(intervaldata,1);
 %  klthreshold=-2.5;
   


   scanIndex=0;
   klvalue=4;
   while klvalue>klthreshold
       scanIndex=scanIndex+1;
       klvalue=KL_calculate(intervaldata(scanIndex,:),template);
   end    
   scan_start=scanIndex;
   Scan_start_afterfilter=Scan_start+scanIndex-1;
   
   scanIndex=IntervalSize+1;
   klvalue=4;
   while klvalue>klthreshold 
       scanIndex=scanIndex-1;
       klvalue=KL_calculate(intervaldata(scanIndex,:),template);
   end   
   scan_end=scanIndex;
   Scan_end_afterfilter=Scan_start+scanIndex-1;
   
   if scan_end>scan_start
      intervalsdata01v2=intervaldata(scan_start:scan_end,:);
      
   else   
      intervalsdata01v2=[];
       Scan_end_afterfilter= Scan_start_afterfilter;
   end   
      
   
   