/*INT and ROUND functions*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/score_data" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
average_score = MEAN (score1, score2, score3);
run;

data scoredata2;
set scoredata1;
average_INT = INT (average_score);
average_ROUND1 = round (average_score, .1);
average_ROUND2 = round (average_score, .01);
run;

proc print DATA = scoredata2;
run;