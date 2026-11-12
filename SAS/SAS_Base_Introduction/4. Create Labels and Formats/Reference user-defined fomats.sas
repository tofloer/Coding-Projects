
/*create data scoredata1*/
proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;  

/*reference/use user-defined formats stored in the Myfmts library*/
libname myfmts "/folders/myfolders/formats" ;
options fmtsearch=(myfmts work library);

proc print data = scoredata1;
	format gender $genderf. averagescore asgroup.;
run;
