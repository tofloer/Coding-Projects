/*Many Obs/subject To One obs/subject using proc transpose*/

libname stuscore "/folders/myfolders";
proc sort data=stuscore.manyper out=manyper;
   by stu_id name;
run;

proc transpose data=manyper 
               out=oneper0;
   by stu_id name;
   id stime ;
   var all_score;/*var to be transposed*/
run;

proc transpose data=manyper prefix=score
               out=oneper1 (drop=_NAME_);
   by stu_id name;
   id stime ;
   /*ID statement names the variable whose formatted values will become 
   the new variable names. */
   var all_score;
run;

/*The PREFIX= option combines the prefix value (score) with the values 
of Stime (1, 2, and 3) to create the three variables score1, score2 and score3. 
PROC TRANSPOSE knew to use the values of Stime to create these variable 
names because it was identified as an ID variable in the procedure.*/
