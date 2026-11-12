/*Proc Means*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;   

PROC MEANS data = scoredata1 MAXDEC = 1 n mean max min;                                   
   VAR score1 score2 score3 averagescore;      
   TITLE 'Summary of Scores';
RUN;

/*By statement*/
PROC SORT DATA = scoredata1 OUT = scoredata2 ;
by Gender ;
run;

PROC MEANS data = scoredata2 MAXDEC = 1 n mean max min;      
   BY gender;                             
   VAR score1 score2 score3 averagescore;      
   TITLE 'Summary of Scores by Gender Using By statement';
RUN;

/*Class statement*/
PROC MEANS data = scoredata1 MAXDEC = 1 n mean max min;      
   class gender;                             
   VAR score1 score2 score3 averagescore;      
   TITLE 'Summary of Scores by Gender Using Class statement';
RUN;

PROC MEANS data = scoredata1 MAXDEC = 1 missing n mean max min;      
   class gender;                             
   VAR score1 score2 score3 averagescore;      
   TITLE 'Summary of Scores by Gender Using Class statement with Missing option';
RUN;

