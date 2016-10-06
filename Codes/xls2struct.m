% this function convert xls to data structrue in matlab

function structure=xls2struct()
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
