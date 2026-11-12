/* Coding Exercise for Section: work with data

coding exercise:
1. Create data set ‘scoredata0’ from importing the excel data ‘score_data_miss777’ 
2. In data step,
a. create a new data set ‘scoredata1’ from ‘scoredata0’ using set statement
b. convert 777 in variables Score1,2,3 to missing (periods) using if… then … else statement
c. calculate average of scores using SAS function MEAN
d. using the claculated averagescore to create Grade categories A, B, C, D, F
3. Print data ‘scoredata1’ using Proc Print
*/

proc import datafile = "/folders/myfolders/score_data_miss777" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;

If score1 = 777 then  score1 = .;
If score2 = 777 then  score2 = .;
If score3 = 777 then  score3 = .;

averagescore = mean (score1, score2, score3);

IF AverageScore >= 90 THEN grade = 'A'; 
Else If averageScore >= 80 then grade = 'B';
Else If averageScore >= 70 then grade = 'C';
Else If averageScore >= 60 then grade = 'D';
Else If 0 =< averageScore < 60 then grade = 'F';
Else grade = ' ';
run;   

PROC PRINT DATA = scoredata1;
RUN;