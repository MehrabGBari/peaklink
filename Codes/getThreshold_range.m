%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     <one line to give the program's name and a brief idea of what it does.>   
%     Copyright (C) <2011>  <Jian Cui>
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


function Threshold=getThreshold_range(XIC,instrument_h_int, de_range,Times_noise_std)
% We determin the threshold for the XIC automatically.
% First we determine noise variance.
XIC=XIC(XIC>0);
if isempty(XIC)==0
dnXIC=wden(XIC,'sqtwolog', 'h','mln',3,'db3');
mvavg=conv(XIC, ones(100,1)/100);
% mvavg=mvavg1;
% cutout values greater then 500
%
mvavg=mvavg(mvavg<100*instrument_h_int/de_range & mvavg>instrument_h_int/de_range);
[histmvavg, intx]=hist(mvavg,100);
[maxhist maxhistid]=max(histmvavg);
meanest=intx(maxhistid);
noisevar=var(XIC-dnXIC);
% noisevar=var(XIC-mvavg1(15:end-15));
noisestd=sqrt(noisevar);

%maxint=max(XIC);

%if maxint/9 > 6*noisestd
%    Threshold=maxint/9;
%else
    Threshold=Times_noise_std*noisestd+meanest;
% Threshold=2*meanest;
else
    Threshold=0;
end   
%end    