/* SAS functions missings  */

proc import datafile = "/folders/myfolders/score_data" 
DBMS = xlsx out = scoredata0 replace ;
run;

DATA scoredata2;
   set scoredata0; 
                                  
   TotalScore = score1 + score2 + score3; /*missing values*/
   AverageScore = TotalScore/3;  /*missing values*/    
   
/*use SAS functions to use only non-missing values for the computation */    
	TotalScore_func = sum (score1, score2, score3);
	AverageScore_func = mean (score1, score2, score3);
	
	Gender_func = UPCASE(gender);                              
	
RUN;

PROC PRINT DATA = scoredata2;
RUN;