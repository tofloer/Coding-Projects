/*proc print*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;   

PROC SORT DATA = scoredata1 OUT = scoredata2 ;
by Gender ;
run;

PROC PRINT DATA = scoredata2;                     
   BY gender;                                 
   SUM AverageScore;                               
   VAR name score1 score2 score3  AverageScore;  
   TITLE 'Proc Print';
RUN;

