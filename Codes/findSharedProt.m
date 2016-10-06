function position = findSharedProt(uniProtMaxquant, pepIntersect)


flagMaxquant = zeros(length(uniProtMaxquant), 1);
for i = 1:length(pepIntersect)
    pos = find(strcmp(uniProtMaxquant, pepIntersect{i})==1);
    if ~isempty(pos)
        flagMaxquant(pos) = 1;
    end
end

position = find(flagMaxquant==1);