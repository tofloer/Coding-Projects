/*Proc Sort*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;   

PROC SORT DATA = scoredata1 OUT = scoredata2 ;
by Gender descending averagescore;
run;

proc print data = scoredata2;
run;

/* options: NODUPKEY  DUPOUT =   */
PROC SORT DATA = scoredata1 OUT = scoredata3 NODUPKEY DUPOUT = extraobs;
by Gender ;
run;

proc print data = scoredata3;
title "sort by gender with NODUPKEY option";
run;
proc print data = extraobs;
title "sort by gender with DUPOUT = extraobs";
run;