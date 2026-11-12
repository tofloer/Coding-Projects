/*Array: creat new vars */

proc import datafile = "/folders/myfolders/score_data_miss999" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1 (drop = i);
set scoredata0;
   ARRAY score_var (3) score1 score2 score3;    
   DO i = 1 TO 3;                       
      IF score_var(i) = 999 THEN score_var(i) =.;   
   END;   
ave_score = mean (score1, score2, score3);
run;

data scoredata2 (drop = i);
set scoredata1;
   array score(3) score1 score2 score3; /*existing vars*/
   array score_diff(3); /*new variables*/
   do i=1 to 3;
      score_diff{i}=score{i} - ave_score;
   end;
run;
