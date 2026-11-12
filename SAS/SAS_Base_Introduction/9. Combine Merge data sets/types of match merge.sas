/*types of match merge*/

/*one to one merge*/

proc import datafile = "/folders/myfolders/score_data_id_gender_only" 
DBMS = xlsx out = scoredata_g replace ;
run;

proc import datafile = "/folders/myfolders/score_data_id_no_gender" 
DBMS = xlsx out = scoredata replace ;
run;

proc sort data = scoredata_g;
by stu_id;
run;

proc sort data = scoredata;
by stu_id;
run;

data OnetoOne;
merge scoredata_g scoredata;
by stu_id;
run;

proc print data = OnetoOne;
title 'one-to-one match merging';
run;

proc print data = scoredata_g;
title 'score data with only gender info.';
run;

proc print data = scoredata;
title 'score data without gender info.';
run;

/*one to many merge*/

proc import datafile = "/folders/myfolders/score_data_id_gender_only_dup" 
DBMS = xlsx out = scoredata_g_dup replace ;
run;

proc import datafile = "/folders/myfolders/score_data_id_no_gender" 
DBMS = xlsx out = scoredata replace ;
run;

proc sort data = scoredata_g_dup;
by stu_id;
run;

proc sort data = scoredata;
by stu_id;
run;

data OnetoMany;
merge scoredata_g_dup scoredata;
by stu_id;
run;

proc print data = OnetoMany;
title 'one-to-many match merging';
run;

proc print data = scoredata_g_dup;
title 'score data with only gender info. with dup.';
run;

proc print data = scoredata;
title 'score data without gender info.';
run;