function plot_matirx(intervalmatrix,number)

colorarray=['r', 'k', 'g', 'b', 'm', 'y','c','b:'];
figure
for kkk=1:size(intervalmatrix,2)
    for i=1:number
        plot(intervalmatrix(:,kkk),colorarray(kkk))
        hold on
    end
    grid on
end



