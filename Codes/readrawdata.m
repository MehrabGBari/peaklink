function [retentiontl1,QCdatal1,peakl,retentiont,MZInt_l1l2]=readrawdata(filename)


% filename='20090608__Orbi6_TaGe_SA_TUMOR_5mix1_02.mzXML';
% fprintf('Please choose the mzXML file')
% [fname,pname]=uigetfile('*.*','open');

% if fname~=0
%     filename=strcat(pname,fname);
%     QCdata=mzxmlread(filename);
% else
%     return;
% end

QCdata=mzxmlread(filename);

totalscan=length(QCdata.scan);
l1count=1;
l2count=1;
for scani=1:totalscan
    if QCdata.scan(scani).msLevel==1
        QCdatal1.scan(l1count)=QCdata.scan(scani);
        l1count=l1count+1;
    else
        QCdatal2.scan(l2count)=QCdata.scan(scani);
        l2count=l2count+1;
    end 
end

totalscan=length(QCdata.scan);
retentiont=zeros(totalscan,2);
for scani=1:totalscan
        if isempty(QCdata.scan(scani).retentionTime(3:end-1))==0
                retentiont(scani,1)=str2num(QCdata.scan(scani).retentionTime(3:end-1));
                retentiont(scani,2)=QCdata.scan(scani).msLevel;
                if isempty(QCdata.scan(scani).precursorMz.value)==0
                    MZInt_l1l2(scani,1)=QCdata.scan(scani).precursorMz.value;
                    MZInt_l1l2(scani,2)=QCdata.scan(scani).precursorMz.precursorIntensity;
                else
                     MZInt_l1l2(scani,1)=0;
                     MZInt_l1l2(scani,2)=0;
                end
        end
end
% save peakvariables03 totalscan retentiont

% save 03QCdata03levelone QCdata03l1 -v7.3;
% save 03QCdata03leveltwo QCdata03l2 -v7.3;
num=0;
for scani=1:length(QCdatal1.scan)
   mz = QCdatal1.scan(scani).peaks.mz(1:2:end);
   int = QCdatal1.scan(scani).peaks.mz(2:2:end);
   
   ID_special=find(mz(1:end-1)-mz(2:end)>=0);
   if ~isempty(ID_special)
       num=num+1;
       mz(ID_special)=[];
       int(ID_special)=[];
   end
   
%    peak=mspeaks(mz,int,'BASE',8);
%    peak=mspeaks(mz,int,'LEVELS',10);
%    peak=mspeaks(mz,int,'NOISEESTIMATOR',NE);
%    peak=mspeaks(mz,int,'MULTIPLIER',1);
%    peak=mspeaks(mz,int,'DENOISING',false);
%    peak=mspeaks(mz,int,'DENOISING',true);
%    peak=mspeaks(mz,int,'PEAKLOCATION',1);
%    peak=mspeaks(mz,int,'FWHHFILTER',0);
%    peak=mspeaks(mz,int,'OVERSEGMENTATIONFILTER',0);
%    peak=mspeaks(mz,int,'HEIGHTFILTER',0);
%    peak=mspeaks(mz,int,'SHOWPLOT',true);
%    peak=mspeaks(mz,int,'SHOWPLOT',false);
   
%    peakl.scan(scani).peak=peak;
    peakl.scan(scani).peak=[mz,int];
end

% scani=1000
%    mz = QCdatal1.scan(scani).peaks.mz(1:2:end);
%    int = QCdatal1.scan(scani).peaks.mz(2:2:end);
% figure
% stem(mz,int)

clear QCdata03;

totalscanl1=length(QCdatal1.scan);
retentiontl1=zeros(totalscanl1,1);
for scani=1:totalscanl1
        if isempty(QCdatal1.scan(scani).retentionTime(3:end-1))==0
                retentiontl1(scani)=str2num(QCdatal1.scan(scani).retentionTime(3:end-1));
        end
end




