/*Array */

proc import datafile = "/folders/myfolders/score_data_miss999" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
   ARRAY score_var (3) score1 score2 score3;    
   DO i = 1 TO 3;                       
      IF score_var(i) = 999 THEN score_var(i) =.;   
   END;       

/*The first time through the DO loop, the variable i has a value of 1 
and the IF statement would read
IF score_var(1) = 999 THEN score_var(1) =.;   
Is the same as 
If score1 = 999 then score1= .;

The second time through, i has a value of 2 
and the IF statement would read 
IF score_var(2) = 999 THEN score_var(2) =.;   
Is the same as 
If score2 = 999 then score2= .;

This continues through all 3 variables in the array
IF score_var(3) = 999 THEN score_var(3) =.;  
Is the same as 
If score3 = 999 then score3= .;
*/                          
RUN;



PROC PRINT data = scoredata1;
title "missing values are shown as periods";
run;

PROC PRINT data = scoredata0;
title "missing values are shown as 999";
run;
