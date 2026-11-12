/*Trim and CATX Function*/

/*create data scoredata0*/
proc import datafile = "/folders/myfolders/chara_data" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
last_name = SCAN(full_name,1, ',: '); /*3 delimiters , : and space*/
first_name = SCAN(full_name,2, ',: '); 
run;

data scoredata2;
set scoredata1;
/*the concatenation operator and Trim function*/
student_name1 = trim(last_name) || ', ' || trim(first_name);
/*the CATX function enables you to concatenate character strings, 
remove leading and trailing blanks, and insert separators*/
length student_name2 $25;
student_name2 = catx(', ',last_name,first_name);
run;



