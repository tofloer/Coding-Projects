/*Substr Function*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/chara_data1" 
DBMS = xlsx out = scoredata0 replace ;
run;


data scoredata1;
set scoredata0;
Exchange=substr(phone,1,3); /*EXTRACT*/
run;

data scoredata2;
set scoredata1;
if exchange='000' then substr(phone,1,3)='408'; /*REPLACE*/
run;

proc print data = scoredata2;
run;
