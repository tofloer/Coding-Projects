/*Index Function*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/chara_data1" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
index_num = index(phone,'408');
run;

data scoredata2;
set scoredata1;
if index_num  = 1 ;
/* if index(phone,'408') = 1;  */
run;


proc print data = scoredata2;
title 'students with phone area code 408';
run;
