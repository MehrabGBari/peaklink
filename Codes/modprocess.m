function [newseq, massdiff, isHeavy]=modprocess(oldseq)
massdiff=0;
newseq=oldseq;
isHeavy=0;

a=strfind(newseq,'C[160.03]');
if ~isempty(a)
massdiff=massdiff+57.0215*numel(a{1});
newseq=strrep(newseq,'C[160.03]','C');
end

a=strfind(newseq,'C[143.00]');
if ~isempty(a)
massdiff=massdiff+(57.0215-17.0265)*numel(a{1});
newseq=strrep(newseq,'C[143.00]','C');
end

a=strfind(newseq,'E[111.03]');
if ~isempty(a)
massdiff=massdiff-18.0106*numel(a{1});
newseq=strrep(newseq,'E[111.03]','E');
end

a=strfind(newseq,'M[147.04]');
if ~isempty(a)
massdiff=massdiff+15.9949*numel(a{1});
newseq=strrep(newseq,'M[147.04]','M');
end

a=strfind(newseq,'Q[111.03]');
if ~isempty(a)
massdiff=massdiff-17.0265*numel(a{1});
newseq=strrep(newseq,'Q[111.03]','Q');
end


a=strfind(newseq,'c[21.01]');
if ~isempty(a)
massdiff=massdiff+4.0085*numel(a{1});
newseq=strrep(newseq,'c[21.01]','');
if numel(a{1})>0
    isHeavy=1;
end
end

[token1, newseq]=strtok(newseq,'.');
newseq=strtok(newseq,'.');
end