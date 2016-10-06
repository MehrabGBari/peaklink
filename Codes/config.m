% config

kl_thresHold = -3.5;


runPreprocessing = 1;
LABEL_TECH = 2; % label = 1 -> SILAC; label = 2 -> super SILAC

PEPTIDE.SEQ = 1;
PEPTIDE.OXIDATION = 44;%50
%PEPTIDE.HLR = 35;
PEPTIDE.CS = 32;
PEPTIDE.MASS = 27;
PEPTIDE.PROTEIN = 28;
PEPTIDE.SCORE = 34;
PEPTIDE.BESTMSMSIDS = 43;%49

MSMS.SEQ = 6;
MSMS.SCANNUM = 3;

OXIDATION.OXIDATIONPROB = 18;
OXIDATION.MZ = 22;
mzXMLPath = 'D:\CoonF4_Mann1_ProfBinMa\MATLAB_Analysis\mzXML\';
% mzXMLMatPath = 'C:\Users\Long\Desktop\Research\SuperSILACAlignment\SuperSalicData\mzXML_mat\';
maxquantPath = 'D:\CoonF4_Mann1_ProfBinMa\MATLAB_Analysis\MaxQuant\';

save_path='..\Preporcessing\';
if ~(exist(save_path,'dir')==7)
    mkdir(save_path);
end
