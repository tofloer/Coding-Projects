/*Array: assign initial values & Creating Temporary Array Elements*/

proc import datafile = "/folders/myfolders/score_data_miss999" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1 (drop= i);
set scoredata0;
   ARRAY score_var (3) score1 score2 score3;    
   DO i = 1 TO 3;                       
      IF score_var(i) = 999 THEN score_var(i) =.;   
   END;   
ave_score = mean (score1, score2, score3);
run;

proc means data = scoredata1;
var score1 score2 score3;    
run; /*79.5, 81.9, 80.8*/

data scoredata2 (drop= i);
set scoredata1;
   array score(3) score1 score2 score3;
   array average{3} (79.5 81.9 80.8); 
   /*assign initial values to average1-3*/
   array score_diff(3);
   do i=1 to 3;
      score_diff{i}=score{i} - average(i);
   end;
run;

proc print data = scoredata2;
run;

/*The variables average1 - 3 are not needed be stored in the data set,
 because they are only for calculating the score_diff1-3. 
 The following program shows you how to create 
 temporary array elements.*/

data scoredata3 (drop= i);
set scoredata1;
   array score(3) score1 score2 score3;
   array average{3} _Temporary_(79.5 81.9 80.80);
   array score_diff(3);
   do i=1 to 3;
      score_diff{i}=score{i} - average(i);
   end;
run;

proc print data = scoredata3;
run;

