function data_information=Saveinformation_1_6data()

fprintf('\n Please choose the testing information xls file: \n')
[fname,pname]=uigetfile('*.*','open');

if fname~=0
    fullname=strcat(pname,fname);
    [a,b,c]=xlsread(fullname);
else
    return;
end
clear a b
[row, col]=size(c);
for i=1:col
    if isnumeric(c{2,i})
        tmp=c(2:end,i);
        ctmp=cell2mat(tmp);
        cc=strcat('structure.',c{1,i},'=ctmp;');
        eval(cc);
    else
        cc=strcat('structure.',c{1,i},'=c(2:end,i);');
        eval(cc);
    end  
end

data_information=structure;
% save data03_information data03_information


% %%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('\n Please choose the input information xls file: \n')
% pepinf=xls2struct();
% 
% % pep01_old=pepinf.peptideseq;
% % chargestate01_old=pepinf.data01ms2chargestate;
% % ms2information01_old=pepinf.data01ms2timepoint;
% % mass01_old=pepinf.data01mass;
% 
% pep01_old=pepinf.peptideseq;
% chargestate01_old=pepinf.data02ms2chargestate;
% ms2information01_old=pepinf.data02ms2timepoint;
% mass01_old=pepinf.data02mass;
% clear pepinf
% 
% for i=1:length(pep01_old)
%     posi01=find(pep01_old{i}=='.');
%     pep01_oldv1{i}=pep01_old{i}(posi01(1)+1:posi01(end)-1);
% end
% for i=1:length(data01_information.peptide)
%     posi01=find(data01_information.peptide{i}=='.');
%     pep01v1{i}=data01_information.peptide{i}(posi01(1)+1:posi01(end)-1);
% end
% Num1=0;
% for i=1:length(pep01_oldv1)
%             TF=strcmp(pep01_oldv1{i},pep01v1);
%             Posi=find(TF==1);
%             if isempty(Posi)==0
% 
%                 [Y,I]=max(data01_information.probability(Posi));
%                 data01_information.retention_time_sec(Posi(I))=ms2information01_old(i);
%                 data01_information.precursor_neutral_mass(Posi(I))=mass01_old(i);
%             end
% end
% data02_information=data01_information;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save data02_information data02_information





