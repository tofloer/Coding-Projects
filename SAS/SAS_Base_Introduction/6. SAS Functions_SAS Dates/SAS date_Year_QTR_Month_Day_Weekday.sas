/*SAS date functions: 
YEAR, QTR, MONTH, and DAY, WEEKDAY functions*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/SAS date" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
sdate = start_date; /*to show the SAS date values for start_date*/
run;

data scoredata2;
set scoredata1;

s_year = year (start_date);
s_qtr = qtr (start_date);
s_month = month (start_date);
s_day = day (start_date);

s_weekday = WEEKDAY(start_date);
run;

data scoredata3;
set scoredata2;
if s_weekday in (1,7); /*indicationg a data error*/
run;
