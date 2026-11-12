/*concatenating*/

proc import datafile = "/folders/myfolders/score_data_id_partial" 
DBMS = xlsx out = scoredata_p replace ;
run;

proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = scoredata replace ;
run;

data CONCAT;
set scoredata_p scoredata;
run;

proc print data = CONCAT;
title 'concatenating';
run;

proc print data = scoredata_p;
title 'partial score data';
run;

proc print data = scoredata;
title 'score data';
run;
