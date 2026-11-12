/* Modify Score dataset with assignment statements */
proc import datafile = "/folders/myfolders/score_data" 
DBMS = xlsx out = scoredata0 replace ;
run;

DATA scoredata1;
   set scoredata0; /*set statement*/
   Grade = 7;  /*new var with a numeric constant */                              
   ScoreType = 'Math Scores';   /*new var with a character constant*/                                           
   TotalScore = score1 + score2 + score3; /*missing values*/
   AverageScore = TotalScore/3;  /*missing values*/    
/*I will show you how to use SAS functions to use only 
non-missing values for the computation in the future tutorials. */    
RUN;
PROC PRINT DATA = scoredata1;
RUN;




