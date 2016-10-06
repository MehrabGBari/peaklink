%%%%% build time model and score seperately from the AR AKL 
%%%%% The warping parameter are A->B A->C and B->C
%%%%% Based on the MS2scan_number_vector, we should get two pairs F1<->B1
%%%%% and F2<->B2.

function TestModel = buildATARAKLModel(Information_Matrix_Total, INTERVAL_Matrix, Parameters, IntervalListBig, retentiontl1_A, retentiontl1_B, retentiontl1_C, retentiont_A, retentiont_B, retentiont_C, XICs_A, XICs_B, XICs_C, PP_Parameter)

disp('Building Models........');
n = size(Information_Matrix_Total,1);
numOfRep = 3;
Total_test_Information_Matrix = cell(n, numOfRep);
TestModel = cell(1,n);



% extract the information from different Information_Matrix_Total,
% INTERVAL_Matrix
for i=1:n
    IsoList=Information_Matrix_Total{i,8}; % get iso distribution
    ID_notzeros=[Information_Matrix_Total{i,4},Information_Matrix_Total{i,12},Information_Matrix_Total{i,20}]~=0; % find the idnetified id
    MS2scan_number_vector=zeros(1,3);
    MS2scan_number_vector(ID_notzeros)=1;    
    
    
    % 
     for filecount=1:3
        if filecount==1 
           TOF_XICs=XICs_A;
           Retent=retentiontl1_A;
           Retentms2=retentiont_A;
       elseif  filecount==2
            TOF_XICs=XICs_B;
            Retent=retentiontl1_B;
            Retentms2=retentiont_B;
       elseif filecount==3
            TOF_XICs=XICs_C;
            Retent=retentiontl1_C;
            Retentms2=retentiont_C;
       end
        
       Interval_after7_combine=IntervalListBig{i,filecount}.intervallist_after7_combine; % get interval list from IntervalListBig
       MS2_ID=INTERVAL_Matrix(i,(filecount-1)*4+2);
       MS2_SC=MS2scan_number_vector(filecount);
       if sum(sum(Interval_after7_combine))==0
           IntervalElutionMatrix{1}=0;
           Interval_after7_combine_Time=0;
           Scan_max=0;
           Time_max_point=0;
       else
           for Interval_Row=1:size(Interval_after7_combine,1)
               % find the elution profile of all possible peaks.
               IntervalElutionMatrix{Interval_Row}=TOF_XICs(Interval_after7_combine(Interval_Row,1):Interval_after7_combine(Interval_Row,2),(i-1)*8+1:i*8);
               Interval_after7_combine_Time{Interval_Row}=Retent(Interval_after7_combine(Interval_Row,1):Interval_after7_combine(Interval_Row,2));
               [a,b]=max(IntervalElutionMatrix{Interval_Row});
               [c,d]=max(a);
               Scan_max(Interval_Row)=b(d);  % the peak position.
               Time_max_point(Interval_Row)=Interval_after7_combine_Time{Interval_Row}(b(d)); % the peak time point.
           end       
       end
        
        Total_test_Information_Matrix{i,filecount}.iso=IsoList;
        Total_test_Information_Matrix{i,filecount}.Interval_after7_combine=Interval_after7_combine;
        Total_test_Information_Matrix{i,filecount}.MS2_ID=MS2_ID;
        Total_test_Information_Matrix{i,filecount}.Interval_after7_combine_Time=Interval_after7_combine_Time;
        Total_test_Information_Matrix{i,filecount}.Scan_max=Scan_max;
        Total_test_Information_Matrix{i,filecount}.Time_max_point=Time_max_point;
        if MS2_SC==0
             Total_test_Information_Matrix{i,filecount}.MS2timePoint=0;
        else
            Total_test_Information_Matrix{i,filecount}.MS2timePoint=Retentms2(MS2_SC);
        end
        Total_test_Information_Matrix{i,filecount}.IntervalElutionMatrix=IntervalElutionMatrix;
        clear Interval_after7_combine_Time IntervalElutionMatrix Interval_after7_combine Scan_max Time_max_point
     end
end


% 
for i=1:size(Information_Matrix_Total,1)
    % J_Val=[INTERVAL_Matrix(i,1),INTERVAL_Matrix(i,5),INTERVAL_Matrix(i,9)];
	J_Val = [1,0,0];
    J_Vec=J_Val~=0;
    J_V=J_Vec*[1;2;4];
    switch J_V
        case 0 %%[A B C]=[0 0 0]

            Row01=0;            Col01=0;            Row02=0;            Col02=0;
            RowElutionProf01{1}=0;            ColElutionProf01{1}=0;
            RowElutionProf02{1}=0;            ColElutionProf02{1}=0;
            
        case 1 %%[A B C]=[1 0 0] (Time warping PP(A)-B and PP(A)-C)
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 2 %%[A B C]=[0 1 0] (Time warping B-PP(A) and PP(B)-C)
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;  
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 3 %%[A B C]=[1 1 0] (Time warping PP(A)-C and PP(B)-C)
            Row01=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,3}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 4 %%[A B C]=[0 0 1] (Time warping C-PP(A) and C-PP(B))
            Row01=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,3}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point; 
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 5 %%[A B C]=[1 0 1] (Time warping PP(A)-B and C-PP(B))
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_BC,Total_test_Information_Matrix{i,2}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point; 
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 6 %%[A B C]=[0 1 1] (Time warping B-PP(A) and C-PP(A))
            Row01=polyval(PP_Parameter.PP_AB,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col01=Total_test_Information_Matrix{i,2}.Time_max_point;
            Row02=polyval(PP_Parameter.PP_AC,Total_test_Information_Matrix{i,1}.Time_max_point);
            Col02=Total_test_Information_Matrix{i,3}.Time_max_point;
            RowElutionProf01=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf01=Total_test_Information_Matrix{i,2}.IntervalElutionMatrix;
            RowElutionProf02=Total_test_Information_Matrix{i,1}.IntervalElutionMatrix;
            ColElutionProf02=Total_test_Information_Matrix{i,3}.IntervalElutionMatrix;
        case 7 %%[A B C]=[1 1 1]    
            %%%%%%%%%%%%%% it is not necessary to do alignment in this
            %%%%%%%%%%%%%% situation
            Row01=0;            Col01=0;            Row02=0;            Col02=0;
            RowElutionProf01{1}=0;            ColElutionProf01{1}=0;
            RowElutionProf02{1}=0;            ColElutionProf02{1}=0;
    end
    TestModel{i}.J_Vec=J_Vec;
    TestModel{i}.Trow01=Row01;
    TestModel{i}.Trow02=Row02;
    TestModel{i}.Tcol01=Col01;
    TestModel{i}.Tcol02=Col02;
    TestModel{i}.Tmatrix01=Row01'*ones(1,length(Col01))-ones(length(Row01),1)*Col01;
    TestModel{i}.Tmatrix02=Row02'*ones(1,length(Col02))-ones(length(Row02),1)*Col02;
    
    [TestModel{i}.ARmatrix01,TestModel{i}.AKLmatrix01]=CalARAKLMatrix(RowElutionProf01,ColElutionProf01);
    [TestModel{i}.ARmatrix02,TestModel{i}.AKLmatrix02]=CalARAKLMatrix(RowElutionProf02,ColElutionProf02);
    for kk=1:3
        TestModel{i}.IntervalListInformation{kk}=IntervalListBig{i,kk};
        TestModel{i}.ImportantInf{kk}=Total_test_Information_Matrix{i,kk};        
    end
    
    clear RowElutionProf01 ColElutionProf01 RowElutionProf02 ColElutionProf02
end
    
for i=1:length(TestModel)
    J_Vec=TestModel{i}.J_Vec;
    J_V=J_Vec*[1;2;4];
    switch J_V
        case 0 %%[A B C]=[0 0 0]
        %%%%%%%%%%%% not the right situation that can align without ms2
        %%%%%%%%%%%% information                    
        case 1 %%[A B C]=[1 0 0] (Time warping PP(A)-B and PP(A)-C)
            Para01=Parameters{1};Para02=Parameters{2};
        case 2 %%[A B C]=[0 1 0] (Time warping B-PP(A) and PP(B)-C)
            Para01=Parameters{1};Para02=Parameters{3};
        case 3 %%[A B C]=[1 1 0] (Time warping PP(A)-C and PP(B)-C)
            Para01=Parameters{2};Para02=Parameters{3};
        case 4 %%[A B C]=[0 0 1] (Time warping C-PP(A) and C-PP(B))
            Para01=Parameters{2};Para02=Parameters{3};
        case 5 %%[A B C]=[1 0 1] (Time warping PP(A)-B and C-PP(B))
            Para01=Parameters{1};Para02=Parameters{3};
        case 6 %%[A B C]=[0 1 1] (Time warping B-PP(A) and C-PP(A))
            Para01=Parameters{1};Para02=Parameters{2};
        case 7 %%[A B C]=[1 1 1]    
            %%%%%%%%%%%%%% it is not necessary to do alignment in this
            %%%%%%%%%%%%%% situation
    end
    if J_V~=0 && J_V~=7
        i
        [TestModel{i}.ATScoreMatrixCorr01,TestModel{i}.ATScoreMatrixNCorr01]=GenATScore(TestModel{i}.Tmatrix01,Para01);
        [TestModel{i}.ATScoreMatrixCorr02,TestModel{i}.ATScoreMatrixNCorr02]=GenATScore(TestModel{i}.Tmatrix02,Para02);
        [TestModel{i}.ARScoreMatrixCorr01,TestModel{i}.ARScoreMatrixNCorr01]=GenARScore(TestModel{i}.ARmatrix01,Para01);
        [TestModel{i}.ARScoreMatrixCorr02,TestModel{i}.ARScoreMatrixNCorr02]=GenARScore(TestModel{i}.ARmatrix02,Para02);
        [TestModel{i}.AKLScoreMatrixCorr01,TestModel{i}.AKLScoreMatrixNCorr01]=GenAKLScore(TestModel{i}.AKLmatrix01,Para01);
        [TestModel{i}.AKLScoreMatrixCorr02,TestModel{i}.AKLScoreMatrixNCorr02]=GenAKLScore(TestModel{i}.AKLmatrix02,Para02);        
        
    end
end