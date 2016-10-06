% Version 1.1 Created by Jianqiu Zhang Feb. 2010
% This function calculates molecular formula, isotopepattern and weight of
% a peptide sequence.
% This function has two input parameters
%       aminosequence: peptide sequence 
%       n:             number of isotopes to be considered.
%It returns these parameters
%       peptideformular: an array with length 5, containing counts of C H N
%                        O and S in the peptide.
%       isotopepattern:   The theoretically calculated isotope pattern with
%                        length n.
%       weight:         The monotopic molecular weight of the peptide.
function [peptideformula,isotopepattern,weight]=aminocalculation_mod(aminosequence,n, modifications)
aminolist=['ILKMFTWVRHANDCEQGPSY'];
%aminoweight=[131.1736 131.1736 146.1882 149.2124 165.1900 119.1197...
%  204.2262 117.1469 174.2017 155.1552 89.0935 132.1184 133.1032 121.1590 147.1299...
%  146.1451 75.0669 115.1310 105.0930 181.1894];
iso1ratio=[1.1078 0.015574 0.3663 0.0372 0.750]/100;
iso2ratio=[0       0        0     0.2004 4.215]/100;
iso3ratio=[0       0        0     0       0.017]/100;
%atomcount=[4.9384 7.7583 1.3577 1.4773 0.0417];   
%averageamino=111.1254;
%atomavgweight=[12.0107 1.00794 14.0067 15.9994 32.065];
atomweight=[12.000  1.0078250321 14.0030740052 15.9949146221 31.97207069];
% The C H N O S
aminoformula=[6 13 1 2 0;
    6 13 1 2 0;
    6 14 2 2 0;
    5 11 1 2 1;
    9 11 1 2 0;
    4 9 1 3 0;
    11 12 2 2 0;
    5 11 1 2 0;
    6 14 4 2 0;
    6 9 3 2 0;
    3 7 1 2 0;
    4 8 2 3 0;
    4 7 1 4 0
    3 7 1 2 1;
    5 9 1 4 0;
    5 10 2 3 0;
    2 5 1 2 0;
    5 9 1 2 0;
    3 7 1 3 0;
    9 11 1 3 0];

al=length(aminosequence);
peptideformula=[0 0 0 0 0];
isotopepattern=zeros(n,1);
weight=0;
for i=1:al
    id=find(aminolist==aminosequence(i));
    if isempty(id)==0
        peptideformula=peptideformula+aminoformula(id,:);
        %weight=weight+aminoweight(id);
    end
end
peptideformula(2)=peptideformula(2)-(al-1)*2;
peptideformula(4)=peptideformula(4)-(al-1);

peptideformula=peptideformula+modifications; %%%%%%%%ADD BY LONG

poisson=peptideformula.*(iso1ratio);
poisson2=peptideformula.*(iso2ratio);

carbon=poisspdf(0:(n-1),poisson(1));
Hydro=poisspdf(0:(n-1),poisson(2));
Nitro=poisspdf(0:(n-1),poisson(3));
Oxy=poisspdf(0:(n-1),poisson(4));
Sul=poisspdf(0:(n-1),poisson(5));

%carbon2=poisspdf(0:(n-1),poisson2(1));
%Hydro2=poisspdf(0:(n-1),poisson2(2));
%Nitro2=poisspdf(0:(n-1),poisson2(3));
Oxy2=poisspdf(0:(n-1),poisson2(4));
Sul2=poisspdf(0:(n-1),poisson2(5));
isotopepattern=conv(carbon,Hydro);
isotopepattern=conv(isotopepattern,Nitro);
isotopepattern=conv(isotopepattern,Oxy);
isotopepattern=conv(isotopepattern,Sul);
isotopepattern=conv(isotopepattern,Oxy2);
isotopepattern=conv(isotopepattern,Sul2);
isotopepattern=isotopepattern(1:n);

weight=sum(peptideformula.*atomweight);
%isotopepattern=isotopepattern/isotopepattern(1);

%I Isoleucine C6H13N1O2 131.1736 CH3-CH2-CH(CH3)-CH(NH2)-COOH
%L Leucine C6H13N1O2 131.1736 (CH3)2-CH-CH2-CH(NH2)-COOH 
%K Lysine C6H14N2O2 146.1882 H2N-(CH2)4-CH(NH2)-COOH 
%M Methionine C5H11N1O2S1 149.2124 CH3-S-(CH2)2-CH(NH2)-COOH 
%F Phenylalanine C9H11N1O2 165.1900 Ph-CH2-CH(NH2)-COOH 
%T Threonine C4H9N1O3 119.1197 CH3-CH(OH)-CH(NH2)-COOH 
%W Tryptophan C11H12N2O2 204.2262 Ph-NH-CH=C-CH2-CH(NH2)-COOH 
%V Valine C5H11N1O2 117.1469 (CH3)2-CH-CH(NH2)-COOH 
%R Arginine C6H14N4O2 174.2017 HN=C(NH2)-NH-(CH2)3-CH(NH2)-COOH 
%H Histidine C6H9N3O2 155.1552 NH-CH=N-CH=C-CH2-CH(NH2)-COOH 
%A Alanine C3H7N1O2  89.0935  CH3-CH(NH2)-COOH
%N Asparagine C4H8N2O3 132.1184 H2N-CO-CH2-CH(NH2)-COOH
%D Aspartate C4H7N1O4 133.1032 HOOC-CH2-CH(NH2)-COOH 
%C Cysteine C3H7N1O2S1 121.1590 HS-CH2-CH(NH2)-COOH
%E Glutamate C5H9N1O4 147.1299 HOOC-(CH2)2-CH(NH2)-COOH 
%Q Glutamine C5H10N2O3 146.1451 H2N-CO-(CH2)2-CH(NH2)-COOH 
%G Glycine C2H5N1O2 75.0669 NH2-CH2-COOH 
%P Proline C5H9N1O2 115.1310 NH-(CH2)3-CH-COOH 
%S Serine C3H7N1O3 105.0930 HO-CH2-CH(NH2)-COOH 
%Y Tyrosine C9H11N1O3 181.1894 HO-Ph-CH2-CH(NH2)-COOH 

