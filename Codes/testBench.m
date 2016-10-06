% Testbench

close all
for i = 1:872
    
    if sum(FinalResult{i}.J_Vec)==3
        [row1, col1] = size(FinalResult{i}.ElutionProfile{1});
        [row2, col2] = size(FinalResult{i}.ElutionProfile{2});
        [row3, col3] = size(FinalResult{i}.ElutionProfile{3});
        
%         if sum(FinalResult{i}.J_Vec)==3
%             score = 'none';
%         else
%             score = [FinalResult{i}.ScoreMatrix01{num}; FinalResult{i}.ScoreMatrix02{num}];
%         end
            
        
        if col1==8&&col2==8&&col3==8
            
            lightProfile = [sum(FinalResult{i}.ElutionProfile{1}(:,1:2), 2); sum(FinalResult{i}.ElutionProfile{2}(:,1:2), 2); sum(FinalResult{i}.ElutionProfile{3}(:,1:2), 2)];
            heavyProfile = [sum(FinalResult{i}.ElutionProfile{1}(:,5:6), 2); sum(FinalResult{i}.ElutionProfile{2}(:,5:6), 2); sum(FinalResult{i}.ElutionProfile{3}(:,5:6), 2)];
            figure
            plot(heavyProfile)
            hold on
            plot(lightProfile, 'r')
            xlabel([num2str(FinalResult{i}.J_Vec), '   |   ', num2str(HLR(i, :))]);
            i
        end
    end
    FinalResult{i}.J_Vec;
end







n = length(proteinSel);
razorProtein = cell(n, 1);
for i = 1:n
    pos = strfind(proteinSel{i}, ';');
    
    if(isempty(pos))
        razorProtein{i} = proteinSel{i};
    else
        razorProtein{i} = proteinSel{i}(1:pos(1)-1);
    end
end


[index1_r1, index2_r1] = uniquePeptideDetection(razorProtein);

ldrUni= log2(HLR(index1_r1, 1)) - log2(HLR(index2_r1, 1));
position = find(ldrUni<100&ldrUni>-100)
std(log2(HLR(index1_r1(position), 1)) - log2(HLR(index2_r1(position), 1)))




median(HLR(pos, 1))
median(HLR(pos, 2))
median(HLR(pos, 3))

figure
subplot(3,1,1)
hist(HLR(pos, 1), 0:0.1:100)
subplot(3,1,2)
hist(HLR(pos, 1), 0:0.1:100)
subplot(3,1,3)
hist(HLR(pos, 1), 0:0.1:100)


std(log2(HLR(pos, 2)) - log2(HLR(pos, 3)))





ms2path='C:\Users\Long\Desktop\SuperSilac\tumor02\combined\txt\peptides.txt';
ms2Info = readtext(ms2path, '\t');

peptideSequence = ms2Info(2:end, 1);
massMaxquant = cell2mat(ms2Info(2:end, 27));
n = length(peptideSequence);
weight = zeros(n, 1);
for i = 1:n
    [peptideformula,isotopepattern,weight(i)]=aminocalculation(peptideSequence{i},4);
    
end
diffMass = weight - massMaxquant;
pos = find(abs(weight - massMaxquant)>1);
i = 1
weight(pos(i))
diffMass(pos(i))
massMaxquant(pos(i))
peptideSequence(pos(i))





