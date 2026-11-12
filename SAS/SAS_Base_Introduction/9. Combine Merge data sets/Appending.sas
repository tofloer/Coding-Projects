/*PROC APPEND*/

proc import datafile = "/folders/myfolders/score_data_id_partial" 
DBMS = xlsx out = scoredata_p replace ;
run;

proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = scoredata replace ;
run;

PROC APPEND BASE=scoredata_p
  DATA=scoredata FORCE;
RUN;

proc print data = scoredata_p;
title 'append';
run;


