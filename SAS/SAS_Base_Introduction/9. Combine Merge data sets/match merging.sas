
/*Match merging */

proc import datafile = "/folders/myfolders/score_data_id_partial_score4" 
DBMS = xlsx out = scoredata_p replace ;
run;

proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = scoredata replace ;
run;

/*ascending*/
proc sort data = scoredata_p;
by stu_id;
run;

proc sort data = scoredata;
by stu_id;
run;

data MMA;
merge scoredata_p scoredata;
by stu_id;
run;

proc print data = MMA;
title 'match merging ascending';
run;

proc print data = scoredata_p;
title 'partial score data';
run;

proc print data = scoredata;
title 'score data';
run;

/*descending*/
proc sort data = scoredata_p;
by descending stu_id;
run;

proc sort data = scoredata;
by descending stu_id;
run;

data MMD;
merge scoredata_p scoredata;
by descending stu_id;
run;

proc print data = MMD;
title 'match merging descending';
run;

proc print data = scoredata_p;
title 'partial score data';
run;

proc print data = scoredata;
title 'score data';
run;