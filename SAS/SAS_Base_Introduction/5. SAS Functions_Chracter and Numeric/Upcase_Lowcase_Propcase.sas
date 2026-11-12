/*case functions*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/chara_data1" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
full_name_LC = lowcase (full_name);
full_name_UC = upcase (full_name);
full_name_PC1 = propcase (full_name_LC, ' ,:');
run;