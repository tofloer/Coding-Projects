/*Assign SAS built-in formats*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;   

data scoredata2;
set scoredata1;
format AverageScore 5.2; /* assign permanent format*/
run;
proc print data = scoredata2;
title "Permanent format for Averagescore";
run;


proc print data = scoredata2;
format AverageScore 4.1; /* assign temporary format*/
title "Temporary format for Averagescore";
run;