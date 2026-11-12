/*SAS date functions: 
today, date, intck functions*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/SAS date" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
today_date1= today();
today_date2= date();
Years=intck('year',start_date,today());
Months=intck('month',start_date,today());
run;

proc print data = scoredata1;
format today_date1 today_date2 date9.;
run;
