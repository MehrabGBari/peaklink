% This script read in the mzXML file and seperate it into levelone scans
% and leveltwo scans.

function [monoXICs,iso1stXICs]=getpepXICmatrix(QCdatal1, pep, chargestate,mass)

% % filename='20090608__Orbi6_TaGe_SA_TUMOR_5mix1_02.mzXML';
% 
% QCdata=mzxmlread(filename);
% 
% totalscan=length(QCdata.scan);
% l1count=1;
% l2count=1;
% for scani=1:totalscan
%     if QCdata.scan(scani).msLevel==1
%         QCdatal1.scan(l1count)=QCdata.scan(scani);
%         l1count=l1count+1;
%     else
%         QCdatal2.scan(l2count)=QCdata.scan(scani);
%         l2count=l2count+1;
%     end 
% end
% 
% totalscan=length(QCdata.scan);
% retentiont=zeros(totalscan,2);
% for scani=1:totalscan
%         if isempty(QCdata.scan(scani).retentionTime(3:end-1))==0
%                 retentiont(scani,1)=str2num(QCdata.scan(scani).retentionTime(3:end-1));
%                 retentiont(scani,2)=QCdata.scan(scani).msLevel;
%         end
% end
% % save peakvariables03 totalscan retentiont
% 
% % save 03QCdata03levelone QCdata03l1 -v7.3;
% % save 03QCdata03leveltwo QCdata03l2 -v7.3;
% clear QCdata03;
% 
% totalscanl1=length(QCdatal1.scan);
% retentiontl1=zeros(totalscanl1,1);
% for scani=1:totalscanl1
%         if isempty(QCdatal1.scan(scani).retentionTime(3:end-1))==0
%                 retentiontl1(scani)=str2num(QCdatal1.scan(scani).retentionTime(3:end-1));
%         end
% end



% supersilac2=mzxmlread(filename);
% totalscan=length(supersilac2.scan);
% l1count=1;
% l2count=1;
% for scani=1:totalscan
%     if supersilac2.scan(scani).msLevel==1
%         supersilac2l1.scan(l1count)=supersilac2.scan(scani);
%         l1count=l1count+1;
%     else
%         supersilac2l2.scan(l2count)=supersilac2.scan(scani);
%         l2count=l2count+1;
%     end 
% end
% 
% totalscan=length(supersilac2.scan);
% retentiont=zeros(totalscan,2);
% for scani=1:totalscan
%         if isempty(supersilac2.scan(scani).retentionTime(3:end-1))==0
%                 retentiont(scani,1)=str2num(supersilac2.scan(scani).retentionTime(3:end-1));
%                 retentiont(scani,2)=supersilac2.scan(scani).msLevel;
%         end
% end
% clear supersilac2;
% 
% %%%%%%%%%%% now extract some peak variables.
% totalscanl1=length(supersilac2l1.scan);
% retentiontl1=zeros(totalscanl1,1);
% for scani=1:totalscanl1
%         if isempty(supersilac2l1.scan(scani).retentionTime(3:end-1))==0
%                 retentiontl1(scani)=str2num(supersilac2l1.scan(scani).retentionTime(3:end-1));
%         end
% end
% next we extract centroid data from levelone scans

%%%%%%%%%%%% get XICs
% n=3;
for i=1:length(pep)
%     aminosequence1=pep{i};
%     [v,p]=find(aminosequence1=='.');
%     aminosequence=aminosequence1(p(1)+1:p(2)-1);    
    pep_chargestate=chargestate(i);
%     [peptideformula,isotopepattern,weight]=aminocalculation(aminosequence,n);
    centermz(i)=(mass(i)+pep_chargestate*1.0073)/pep_chargestate;
    isotope1st_mz(i)=(mass(i)+1.00335+pep_chargestate*1.0073)/pep_chargestate;
    isotope2nd_mz(i)=(mass(i)+2*1.00335+pep_chargestate*1.0073)/pep_chargestate;
end

tolerance=5;  %% 5ppm
monoXICs=getXICs(QCdatal1,centermz,tolerance);
iso1stXICs=getXICs(QCdatal1,isotope1st_mz,tolerance);
% iso2ndXICs=getXICs(QCdatal1,isotope2nd_mz,tolerance);

% monoXICdetails=getrawXICs(QCdatal1,centermz,tolerance);
% iso1stXICdetails=getrawXICs(QCdatal1,isotope1st_mz,tolerance);



















  
