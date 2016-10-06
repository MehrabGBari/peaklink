function [Total_intervallist_new L]=removeShortIntervals(Total_intervallist_new,minlength)
   intervalLength=Total_intervallist_new(:,2)-Total_intervallist_new(:,1);
   idskeep=intervalLength>=minlength;
   if isempty(idskeep)==0
     Total_intervallist_new=Total_intervallist_new(idskeep,:);
     L=sum(idskeep);
   else
     Total_intervallist_new=[0 0];
     L=0;
   end  
       