% logging

dateInfo = date;

fid = fopen([save_path, ExperimentID,'_log.txt'],'wt');


fprintf(fid, '%s\n\n', [ExperimentID,'----',dateInfo, '----']);

fprintf(fid, '%s\t', 'time used:');
fprintf(fid, '%s\n', num2str(timeConsumed));


fprintf(fid, '%s\n', '------------FILE NAMES------------');
fprintf(fid, '\t%s\n', ['1.    ', mzXMLfileName1]);
fprintf(fid, '\t%s\n', ['2.    ', mzXMLfileName2]);
fprintf(fid, '\t%s\n\n', ['3.    ', mzXMLfileName3]);

fprintf(fid, '%s\n', '------------PATHS------------');
fprintf(fid, '\t%s\n', path1);
fprintf(fid, '\t%s\n', path2);
fprintf(fid, '\t%s\n', path3);
fprintf(fid, '\t%s\n', MzxmlFile_path1);
fprintf(fid, '\t%s\n', MzxmlFile_path2);
fprintf(fid, '\t%s\n', MzxmlFile_path3);


fprintf(fid, '%s\n', '------------RESULT----ErrorRate-----');

fprintf(fid, '\t%s\n', '--------1 vs 2---------');

fprintf(fid, '\t\t%s\n', 'num of intetsection: ');
fprintf(fid, '\t\t\t%s\n', num2str(nIntersectErrRate(1)));
fprintf(fid, '\t\t%s\n', 'num of trainingsets: ');
fprintf(fid, '\t\t\t%s\n', num2str(200));
fprintf(fid, '\t\t%s\n', 'scores: AR, AT, AR + AT');
fprintf(fid, '\t\t\t%s\t', num2str(scoreErrRate(1, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(scoreErrRate(1, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(scoreErrRate(1, 3)));
fprintf(fid, '\t\t%s\n', 'num of testing sets: ');
fprintf(fid, '\t\t\t%s\t', num2str(nTestingErrRate(1, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(nTestingErrRate(1, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(nTestingErrRate(1, 3)));


fprintf(fid, '\t%s\n', '--------1 vs 3---------');
fprintf(fid, '\t\t%s\n', 'num of intetsection: ');
fprintf(fid, '\t\t\t%s\n', num2str(nIntersectErrRate(2)));
fprintf(fid, '\t\t%s\n', 'num of trainingsets: ');
fprintf(fid, '\t\t\t%s\n', num2str(200));
fprintf(fid, '\t\t%s\n', 'scores: AR, AT, AR + AT');
fprintf(fid, '\t\t\t%s\t', num2str(scoreErrRate(2, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(scoreErrRate(2, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(scoreErrRate(2, 3)));
fprintf(fid, '\t\t%s\n', 'num of testing sets: ');
fprintf(fid, '\t\t\t%s\t', num2str(nTestingErrRate(2, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(nTestingErrRate(2, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(nTestingErrRate(2, 3)));


fprintf(fid, '\t%s\n', '--------total number of entry: ---------');
fprintf(fid, '\t\t%s\t', 'Data 1:');
fprintf(fid, '%s\n', num2str(nData1));

fprintf(fid, '\t\t%s\t', 'Data 2:');
fprintf(fid, '%s\n', num2str(nData2));

fprintf(fid, '\t\t%s\t', 'Data 3:');
fprintf(fid, '%s\n', num2str(nData3));

fprintf(fid, '\t\t%s\t', 'Union:');
fprintf(fid, '%s\n', num2str(nUnion));


fprintf(fid, '%s\n', '------------RESULT----Detection Rate-----');


fprintf(fid, '\t%s\n', '--------1 vs 2---------');

fprintf(fid, '\t\t%s\n', 'num of intetsection: ');
fprintf(fid, '\t\t\t%s\n', num2str(nIntersectDetRate(1)));
fprintf(fid, '\t\t%s\n', 'num of trainingsets: ');
fprintf(fid, '\t\t\t%s\n', num2str(200));
fprintf(fid, '\t\t%s\n', 'scores: AR, AT, AR + AT');
fprintf(fid, '\t\t\t%s\t', num2str(scoreDetRate(1, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(scoreDetRate(1, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(scoreDetRate(1, 3)));
fprintf(fid, '\t\t%s\n', 'num of testing sets: ');
fprintf(fid, '\t\t\t%s\t', num2str(nTestingDetRate(1, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(nTestingDetRate(1, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(nTestingDetRate(1, 3)));


fprintf(fid, '\t%s\n', '--------1 vs 3---------');
fprintf(fid, '\t\t%s\n', 'num of intetsection: ');
fprintf(fid, '\t\t\t%s\n', num2str(nIntersectDetRate(2)));
fprintf(fid, '\t\t%s\n', 'num of trainingsets: ');
fprintf(fid, '\t\t\t%s\n', num2str(200));
fprintf(fid, '\t\t%s\n', 'scores: AR, AT, AR + AT');
fprintf(fid, '\t\t\t%s\t', num2str(scoreDetRate(2, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(scoreDetRate(2, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(scoreDetRate(2, 3)));
fprintf(fid, '\t\t%s\n', 'num of testing sets: ');
fprintf(fid, '\t\t\t%s\t', num2str(nTestingDetRate(2, 1)));
fprintf(fid, '\t\t\t%s\t', num2str(nTestingDetRate(2, 2)));
fprintf(fid, '\t\t\t%s\n', num2str(nTestingDetRate(2, 3)));


fprintf(fid, '\t%s\n', '--------total number of entry: ---------');
fprintf(fid, '\t\t%s\t', 'Data 1:');
fprintf(fid, '%s\n', num2str(nData1));

fprintf(fid, '\t\t%s\t', 'Data 2:');
fprintf(fid, '%s\n', num2str(nData2));

fprintf(fid, '\t\t%s\t', 'Data 3:');
fprintf(fid, '%s\n', num2str(nData3));

fprintf(fid, '\t\t%s\t', 'Union:');
fprintf(fid, '%s\n', num2str(nUnion));


fclose(fid);
