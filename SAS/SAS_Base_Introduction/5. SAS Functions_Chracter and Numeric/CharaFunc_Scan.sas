/*Scan function*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/chara_data" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
last_name = SCAN(full_name,1, ',: '); /*3 delimiters , : and space*/
run;

proc print data = scoredata1;
run;
