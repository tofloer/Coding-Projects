/*Section: Error handling
Coding Exercise:
data:
score_data_miss777.xlsx ---> output SAS data set sdm0
instruction:
1. using PUTLOG to detect the possible logic error of the following data step program:
DATA sdm1;
   set sdm0;                                           
   AverageScore = mean (score1, score2, score3);  
   If averagescore <60; 
run;
2. correct the error
*/
proc import datafile = "/folders/myfolders/score_data_miss777" 
DBMS = xlsx out = sdm0 replace ;
run;

DATA sdm1;
   set sdm0;                                           
   AverageScore = mean (score1, score2, score3);  
   PUTLOG 'the score infomation: ' Name= Score1= Score2= Score3= AverageScore= 5.2;
   If averagescore <60 ; 
run;

data sdm2;
set sdm0;
ARRAY score_var (3) score1 score2 score3;    
   DO i = 1 TO 3;                       
      IF score_var(i) = 777 THEN score_var(i) =.;   
   END;     
AverageScore = mean (score1, score2, score3);  
   PUTLOG 'the score infomation: ' Name= Score1= Score2= Score3= AverageScore= 5.2;
   If averagescore <60 ;
run;
