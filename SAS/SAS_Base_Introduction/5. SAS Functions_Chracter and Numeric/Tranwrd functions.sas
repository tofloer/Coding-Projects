/*TRANWRD Function*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/chara_data1" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
length phone_update $10;
phone_update = TRANWRD(phone,'000','408');
run;
