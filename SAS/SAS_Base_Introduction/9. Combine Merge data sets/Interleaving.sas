/*Interleaving */

proc import datafile = "/folders/myfolders/score_data_id_partial" 
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

data interleave;
set scoredata_p scoredata;
by stu_id;
run;

proc print data = interleave;
title 'interleaving';
run;

proc print data = scoredata_p;
title 'partial score data';
run;

proc print data = scoredata;
title 'score data';
run;