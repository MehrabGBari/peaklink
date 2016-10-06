%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     <one line to give the program's name and a brief idea of what it does.>   
%     Copyright (C) <2011>  <Jian Cui>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Affero General Public License as
%     published by the Free Software Foundation, either version 3 of the
%     License, or (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Affero General Public License for more details.
% 
%     You should have received a copy of the GNU Affero General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [data01v1,I,P01_once]=filteroverlapp(data01_information)
% data01=data01_information;
[Y,I]=sort(data01_information.probability,'descend');
data01.peptide=data01_information.peptide(I);
data01.probability=data01_information.probability(I);
data01.assumed_charge=data01_information.assumed_charge(I);
data01.calc_neutral_pep_mass=data01_information.calc_neutral_pep_mass(I);
data01.remarks=data01_information.remarks(I);
data01.MZratio=data01_information.MZratio(I);
data01.retention_time_sec=data01_information.retention_time_sec(I);
data01.protein=data01_information.protein(I);
[Pepunique01,P01_once]=unique(data01.peptide,'first');

% for i=1:length(data01.peptide)
%     if data01.index(i)~=0
%         pepseq=data01.peptide{i};
%         TF=strcmp(pepseq,data01.peptide);
%         posi=find(TF==1);
%         if length(posi)>1
%             Prob=data01.probability(posi);
%             [v,p]=max(Prob);
%             posi(p)=[];
%             data01.index(posi)=0;
%         end
%     end
% end
% 
% P01=find(data01.index~=0);
% P01_once=data01.index(P01);

data01v1.peptide=data01.peptide(P01_once);
data01v1.probability=data01.probability(P01_once);
data01v1.assumed_charge=data01.assumed_charge(P01_once);
data01v1.remarks=data01.remarks(P01_once);
data01v1.retention_time_sec=data01.retention_time_sec(P01_once);
data01v1.calc_neutral_pep_mass=data01.calc_neutral_pep_mass(P01_once);
data01v1.MZratio=data01.MZratio(P01_once);
data01v1.protein=data01.protein(P01_once);






