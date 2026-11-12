/*user defined formats*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;   

PROC FORMAT;                            
   VALUE $genderf 	'm' = 'Male'              
                	'f' = 'Female';   
   VALUE asgroup    0-<60 = 'F'
   					60-<70 = 'D'
   					70-<80 = 'C'  
   					80-<90 = 'B'
   					90-High = 'A'
   					Other = 'Missing';
run;
   
PROC PRINT DATA = scoredata1;
   FORMAT gender $genderf. averagescore asgroup.;
   TITLE 'Results Printed with User-Defined Formats';
RUN;

PROC PRINT DATA = scoredata1;
   TITLE 'Results Printed without User-Defined Formats';
RUN;

