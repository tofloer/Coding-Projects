/* subset your data*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

DATA scoredata_all;
   set scoredata0; 
   length take $20; /*code update: adding enougth length to Take*/
   /*take*/
	IF score1 NE . AND score2 NE .  AND score3 NE .  THEN take = 'complete';
	else take = 'incomplete';
	
run;

/* subsetting IF*/
data scoredata_subset_if;
set scoredata_all;
if take = 'complete';
run;
proc print data = scoredata_subset_if;
run;

/*Delete*/
data scoredata_subset_delete;
set scoredata_all;
if take NE 'complete' then delete;
run;
proc print data = scoredata_subset_delete;
run;


    