/*One to one merge*/

proc import datafile = "/folders/myfolders/score_data_id_partial" 
DBMS = xlsx out = scoredata_p replace ;
run;

proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = scoredata replace ;
run;

data onetoone;
set scoredata_p;
set scoredata;
run;

proc print data = onetoone;
title 'one to one merge';
run;

proc print data = scoredata_p;
title 'partial score data';
run;

proc print data = scoredata;
title 'score data';
run;




