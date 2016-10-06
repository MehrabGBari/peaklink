function [Ground_trainingtimepoint01,Ground_trainingtimepoint02,Ground_trainingtimepoint03,...
 PP12,PP13,PP21,PP23,PP31,PP32]=GenerateGwarpingParameter(Training_matrix,...
    groundtruthinterval01_final,groundtruthinterval02_final,groundtruthinterval03_final,...
    IntervalList01_final,IntervalList02_final,IntervalList03_final,...
    monoXICs01_final,monoXICs02_final,monoXICs03_final,...
    retentiont01l1,retentiont02l1,retentiont03l1)
for i=1:size(Training_matrix,1)
    
    id_in01=Training_matrix(i,8);
    id_in02=Training_matrix(i,16);
    id_in03=Training_matrix(i,24);
    
    id_intervalindex01=groundtruthinterval01_final(id_in01,2);
    id_intervalindex02=groundtruthinterval02_final(id_in02,2);
    id_intervalindex03=groundtruthinterval03_final(id_in03,2);
    
    interval01=IntervalList01_final{id_in01}.intervallist(id_intervalindex01,:);
    interval02=IntervalList02_final{id_in02}.intervallist(id_intervalindex02,:);
    interval03=IntervalList03_final{id_in03}.intervallist(id_intervalindex03,:);
    
    Ground_elution_peak01=monoXICs01_final(interval01(1):interval01(2),id_in01);
    Ground_elution_peak02=monoXICs02_final(interval02(1):interval02(2),id_in02);
    Ground_elution_peak03=monoXICs03_final(interval03(1):interval03(2),id_in03);
    
    %%%%%% retentiont format is different
    Ground_trainingtimepoint01(i)=sum(retentiont01l1(interval01))/2;
    Ground_trainingtimepoint02(i)=sum(retentiont02l1(interval02))/2;
    Ground_trainingtimepoint03(i)=sum(retentiont03l1(interval03))/2;
    
    
end

PP12=polyfit(Ground_trainingtimepoint01,Ground_trainingtimepoint02,4);
PP13=polyfit(Ground_trainingtimepoint01,Ground_trainingtimepoint03,4);
PP21=polyfit(Ground_trainingtimepoint02,Ground_trainingtimepoint01,4);
PP23=polyfit(Ground_trainingtimepoint02,Ground_trainingtimepoint03,4);
PP31=polyfit(Ground_trainingtimepoint03,Ground_trainingtimepoint01,4);
PP32=polyfit(Ground_trainingtimepoint03,Ground_trainingtimepoint02,4);
% figure
% plot(Ground_trainingtimepoint01,Ground_trainingtimepoint03,'r.')
% XXX=1:50:14000; YYY=polyval(PP13,XXX);
% hold on; plot(XXX,YYY);