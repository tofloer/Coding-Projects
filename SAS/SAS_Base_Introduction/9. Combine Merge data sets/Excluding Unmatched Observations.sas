
/*Match merging -- Excluding Unmatched Observations*/

proc import datafile = "/folders/myfolders/score_data_id_partial_score4" 
DBMS = xlsx out = scoredata_p replace ;
run;

proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = scoredata replace ;
run;

proc sort data = scoredata_p;
by stu_id;
run;

proc sort data = scoredata;
by stu_id;
run;

data MMA;
merge scoredata_p(in = A) scoredata (in = B);
/*Creating Temporary IN= Variables a and b*/
/*if scoredata_p contributed data to the current observation,
  A = 1, otherwise A = 0, etc*/
by stu_id;
if A = 1 and B = 1;
run;

proc print data = MMA;
title 'Selecting Matching Observations';
run;

proc print data = scoredata_p;
title 'partial score data';
run;

proc print data = scoredata;
title 'score data';
run;
