function index = findUniquePep(proteinList)

dd = 1;
clear index;
for i = 1:length(proteinList)
    pos = strfind(proteinList{i}, ';');
    if isempty(pos)
        index(dd) = i;
        dd = dd + 1;
    end
end