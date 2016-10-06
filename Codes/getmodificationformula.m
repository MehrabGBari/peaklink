function M=getmodificationformula(seq)
M=[0 0 0 0 0];
% C H N O S
a=strfind(seq,'C[160.03]');
M=M+[2 3 1 1 0]*length(a);

a=strfind(seq,'C[143.00]');
M=M+[0 -3 -1 0 0]*length(a);

a=strfind(seq,'E[111.03]');
M=M+[0 -2 0 -1 0]*length(a);

a=strfind(seq,'M[147.04]');
M=M+[0 0 0 1 0]*length(a);

a=strfind(seq,'Q[111.03]');
M=M+[0 -3 -1 0 0]*length(a);

a=strfind(seq,'n[43.02]');
M=M+[2 2 0 1 0]*length(a);
