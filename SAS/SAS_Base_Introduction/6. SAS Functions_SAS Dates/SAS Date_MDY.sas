/*SAS date MDY*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/SAS date" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
Bdate = mdy (month, day, year);
fix_date1 = mdy (8, 31,1950);
fix_date2 = mdy (8, 31,2050);
fix_date3 = mdy (8, 31,50);
run;

/*
The YEARCUTOFF= option specifies the first year of a 100 year 
    window in which all 2-digit years are assumed to be. For 
    example, if the YEARCUTOFF= option is set to 1920, all 2-digit 
    years are assumed to be in the period 1920 through 2019. This 
    means that two-digit years from 20-99 will be assigned a 
    century prefix of '19' and all 2-digit years from 00-19 will 
    have a century prefix of '20'. 
*/

proc print data = scoredata1;
format bdate date9. fix_date1 date9. fix_date2 date9. fix_date3 date9.;
run;
