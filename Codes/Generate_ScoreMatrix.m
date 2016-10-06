function         [ScoreMatrix01, ScoreMatrix02]=Generate_ScoreMatrix(TestModel_spec)


ATMatrix01=log(TestModel_spec.ATScoreMatrixCorr01./TestModel_spec.ATScoreMatrixNCorr01);
% ATMatrix01(TestModel_spec.ATScoreMatrixCorr01<=TestModel_spec.ATScoreMatrixNCorr01) = -100;


ARMatrix01=log(TestModel_spec.ARScoreMatrixCorr01./TestModel_spec.ARScoreMatrixNCorr01);
% ARMatrix01(TestModel_spec.ARScoreMatrixCorr01<=TestModel_spec.ARScoreMatrixNCorr01) = -100;

AKLMatrix01=log(TestModel_spec.AKLScoreMatrixCorr01./TestModel_spec.AKLScoreMatrixNCorr01);
% AKLMatrix01(TestModel_spec.AKLScoreMatrixCorr01<=TestModel_spec.AKLScoreMatrixNCorr01) = -100;


ATMatrix02=log(TestModel_spec.ATScoreMatrixCorr02./TestModel_spec.ATScoreMatrixNCorr02);
% ATMatrix02(TestModel_spec.ATScoreMatrixCorr02<=TestModel_spec.ATScoreMatrixNCorr02) = -100;

ARMatrix02=log(TestModel_spec.ARScoreMatrixCorr02./TestModel_spec.ARScoreMatrixNCorr02);
% ARMatrix02(TestModel_spec.ARScoreMatrixCorr02<=TestModel_spec.ARScoreMatrixNCorr02) = -100;


AKLMatrix02=log(TestModel_spec.AKLScoreMatrixCorr02./TestModel_spec.AKLScoreMatrixNCorr02);
% AKLMatrix02(TestModel_spec.AKLScoreMatrixCorr02<=TestModel_spec.AKLScoreMatrixNCorr02) = -100;


ScoreMatrix01{1}=ATMatrix01;
ScoreMatrix01{2}=ARMatrix01;
ScoreMatrix01{3}=AKLMatrix01;
ScoreMatrix01{4}=ATMatrix01+ARMatrix01;
ScoreMatrix01{5}=ATMatrix01+AKLMatrix01;
ScoreMatrix01{6}=ARMatrix01+AKLMatrix01;
ScoreMatrix01{7}=ATMatrix01+ARMatrix01+AKLMatrix01;


ScoreMatrix02{1}=ATMatrix02;
ScoreMatrix02{2}=ARMatrix02;
ScoreMatrix02{3}=AKLMatrix02;
ScoreMatrix02{4}=ATMatrix02+ARMatrix02;
ScoreMatrix02{5}=ATMatrix02+AKLMatrix02;
ScoreMatrix02{6}=ARMatrix02+AKLMatrix02;
ScoreMatrix02{7}=ATMatrix02+ARMatrix02+AKLMatrix02;