/*section: Combining SAS data sets
coding exercise:

data:
score_data_id_class.xls: student-level info
class_info.xls: class-level info

instruction:
merge class information to student-level data by variable class. 
Only the observations with class-level information will
be kept in the merged data.
*/


proc import datafile = "/folders/myfolders/score_data_id_class" 
DBMS = xlsx out = sd0 replace ;
run;
proc import datafile = "/folders/myfolders/class_info" 
DBMS = xlsx out = c0 replace ;
run;

proc sort data = sd0;
by class;
run;
proc sort data = c0;
by class;
run;

data m0;
merge c0 (in = inC) sd0 (in = inS);
by class;
if inC;
run;