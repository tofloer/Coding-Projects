/*SAS date functions: 
DATDIF, YRDIF functions*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/SAS date" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
year_diff1 = YRDIF(start_date,today(),'30/360');
year_diff2 = YRDIF(start_date,'14sep2015'd,'ACT/ACT');
year_diff3 = YRDIF(start_date,today(),'ACT/360');
year_diff4 = YRDIF(start_date,'14sep2015'd,'ACT/365');

day_diff1 = datDIF(start_date,today(),'30/360');
day_diff2 = datDIF(start_date,'14sep2015'd,'ACT/ACT');

run;

proc print data = scoredata1;
run;