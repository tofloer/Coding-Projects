/* SAS SQL Essentials: Proc SQL

coding exercise:
data: 
allscore.xlsx --> import to create SAS data set allScore0
score_data_id.xlsx --> import to create SAS data set Score_data0

instruction:
1. using SELECT to create a new table 'allScore1' with variables stuid, math from table 'allScore0'
2. using SELECT to create a new table 'Score_data1' with variables stu_id, score2 from table 'Score_data0'
3. join two tables Score_data1 and allScore1 on stuid = stu_id Using keyword JOIN
a. keep only the matching student 
b. keep all student from table allScore1 
c. include all student in two tables*/

/*1*/
proc import datafile = "/folders/myfolders/allscore" 
DBMS = xlsx out = allscore0 replace ;
run;

proc sql;
   create table allscore1 as
   select stuid,
          math
   from allscore0;
quit;

proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = score_data0 ;
run;

/*2*/
proc sql;
   create table score_data1 as
   select stu_id,
          score2
   from score_data0;
quit;

/*3*/
/*a. keep only the matching student*/
proc sql;
   CREATE TABLE inner AS
   SELECT stu_id, stuid, math, score2
   FROM score_data1 AS sd1 INNER JOIN  allscore1 AS as1
   ON stu_id =  stuid;
quit;

/*b. keep all student from table Allscore1 */
proc sql;
   CREATE TABLE from_allscore1 AS
   SELECT stu_id, stuid, math, score2
   FROM score_data1 AS sd1 RIGHT JOIN  allscore1 AS as1
   ON stu_id =  stuid;
quit;

/*c. include all student in two tables*/
proc sql;
   CREATE TABLE two AS
   SELECT stu_id, stuid, math, score2
   FROM score_data1 AS sd1 FULL JOIN  allscore1 AS as1
   ON stu_id =  stuid;
quit;
