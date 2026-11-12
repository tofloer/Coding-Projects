/*Converting a Data Set with Several Observations per Subject 
to a Data Set with One Observation per Subject: 
Using a DATA Step*/

libname stuscore "/folders/myfolders";
proc sort data=stuscore.manyper out=manyper;
   by stu_id stime; /*stime is the number of scores*/
run;

data stuscore.oneper;
   set manyper;
   by stu_id stime;
   array score{3}; /*creating score1, score2 and score3*/
   retain score1-score3;
/*RETAIN statement:::
You need to retain these three variables because they do not 
come from a SAS data set and are, by default, set equal to a missing 
value for each iteration of the DATA step. The RETAIN statement 
prevents this from happening.*/
   if first.stu_id then call missing(of score1-score3);
/* FIRST.variable and LAST.variable:::
If you use a by statement along with a set statement in a data step then 
SAS creates two automatic variables, FIRST.variable and LAST.variable, 
where variable is the name of the by variable.  
FIRST.variable has a value 1 for the first observation in the by group and 0 
for all other observations in the by group. 
LAST.variable has a value 1 for the last observation in the by group and 
0 for all other observations in the by group.
You can take actions conditionally, based on whether you are processing the first 
or the last observation in a BY group.
These temporary variables are available for DATA step programming but are n
ot added to the output data set.*/

/*when you start processing the first visit for each subject/student, you set 
the three values of score to missing. Otherwise, a subject with fewer 
than three scores may wind up with a score from the previous student. */

/*The CALL MISSING routine :::
can set any number of numeric and/or character 
values to missing at one time. */

   score{stime} = all_score;
   if last.stu_id then output;
   keep name stu_id score1-score3;
run;

proc print data = stuscore.manyper;
title 'multiple records per student';
run;
proc print data = stuscore.oneper;
title 'one record per student';
run;

