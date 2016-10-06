%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     <one line to give the program's name and a brief idea of what it does.>   
%     Copyright (C) <2010>  <Jianqiu Zhang>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Affero General Public License as
%     published by the Free Software Foundation, either version 3 of the
%     License, or (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Affero General Public License for more details.
% 
%     You should have received a copy of the GNU Affero General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version 1.1 Created by Jianqiu Zhang Feb. 2010
% This function returns a list a intervals of the XIC above the threshold.
% The intervallist will be sorted according to the maximum intensity value
% during each interval.
% Inputs: 
%     XIC: a column sigle array;
%     th: threshold to be applied to xic;
%     mininterval: Mimimu length of an interval.
%     maxnointervals: Maximum number of intervals considered.
% Outputs:
%     intervallist: a two column matrix. Each row contains the start and the
%     end of the interval. The number of the intervals is  maxinointervals or smaller.
%     If no intervals are detected, it returns an [0 0]. When
%     calling this function, always check if the returned list is empty.
%     For example:
%       intervallist=intervaldetection(xic,th,3,5);
%       if intervallist(1,1)==0
%              .....
%       end

function [intervallist]=TOF_intervaldetectionv1(xic_original,xic,mininterval,maxnointervals,min_interval_length,th)

ll=length(xic);
intervalmap=zeros(ll,1);
if isempty(th)==1
    th=20;
end
ids=find(xic>=th);
%intervallist=zeros(maxnointervals,2);
if isempty(ids)==0
 intervalmap(ids)=1;
 
 % added the following line to correct for when the last point in xic is
 % above the threshold...this phenomenon will crash MATLAB
 
 intervalmap(ll) = 0;

 shiftmap=circshift(intervalmap,1);
 diff=intervalmap-shiftmap;
 intervalstart=find(diff==1);
 diff(intervalstart)=0;
 diff=circshift(diff,-1);
 intervalend=find(diff==-1);
 linterval=length(intervalstart);
 intervalend=intervalend(1:linterval);
 
 intervallength=intervalend-intervalstart+1;
 iids=find(intervallength>mininterval);
 intervalstart=intervalstart(iids);
 intervalend=intervalend(iids);
 linterval=length(intervalstart);
 % added the next lines to correct for having no intervals long enough to
 % pass the mininterval test...this phenomenon will also crash MATLAB
 
 if (isempty(intervalstart) || isempty(intervalend))
     intervallist=[ 0  0];
     return
 end
 heightv=zeros(linterval,1);
 for i=1:linterval
     heightv(i)=max(xic(intervalstart(i):intervalend(i)));
 end    
 
 %intervallength=intervallength(iids);
 %[dummy, sortid]=sort(intervallength,'descend');
 [dummy, sortid]=sort(heightv,'descend');
 intervalstart=intervalstart(sortid);
 intervalend=intervalend(sortid);
 
   if (linterval>maxnointervals)
     intervallist=[intervalstart(1:maxnointervals) intervalend(1:maxnointervals)];   
   else  
     intervallist=[intervalstart intervalend];
   end
 else
    intervallist=[0  0];
 end   


ID_interval_short=find(intervallist(:,2)-intervallist(:,1)<=min_interval_length & intervallist(:,2)-intervallist(:,1)>0);
intervallist(ID_interval_short,:)=[];
if isempty(intervallist)
     intervallist=[0  0];
end

 
 
 
 
