/* conditonal processing using If Then statements*/

proc import datafile = "/folders/myfolders/score_data" 
DBMS = xlsx out = scoredata0 replace ;
run;

DATA scoredata_ifthen;
   set scoredata0; 
   /*using SAS functions*/
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
	/*using If...Then... statements*/
	If gender = 'm' then gender_num = 1;
	IF score1 NE . AND score2 NE .  AND score3 NE .  THEN take = 'complete';
	IF AverageScore >= 90 THEN DO;
       grade = 'A';
       pass = 'pass';
	END;  
run;                            		
	
/* there are missing values in new vars created using IF THEN statements,
it is because the code only addressed one condition in each new variable, 
in the next tutorial I will show you to use IF-THEN/ELSE Statements to 
complete creating new variables*/
                           
	


