/*Converting a Data Set with One Observation per Subject to 
a Data Set with Several Observations per Subject: 
Using a DATA Step
*/
proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = oneper replace ;
run;

libname stuscore "/folders/myfolders";

data stuscore.manyper;
   set oneper;
   array score{3}score1 score2 score3;
/*The score array has three elements: score1 score2 score3*/ 
   do stime = 1 to 3;
      all_score = score{stime}; 
/*when stime = 1, all_score = score1; when stime =2, all_score = score2, etc*/
	  if all_score ne . then output;
   end;
   keep name stu_id all_score stime;
run;

proc print data = oneper;
title 'one obs per student'; 
run;
proc print data = stuscore.manyper;
title 'several obs per student'; 
run;





