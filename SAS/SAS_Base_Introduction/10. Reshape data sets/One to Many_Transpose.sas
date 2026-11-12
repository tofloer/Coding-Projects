/* One to Many using Proc Transpose*/

proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = oneper replace ;
run;

proc sort data = oneper ;
by stu_id;
run;

/*PROC TRANSPOSE takes an input data set and outputs a data set where the 
original rows become columns and the original columns become rows. */
proc transpose data=oneper
               out=manyper;
   by stu_id;
   /*BY statement that performs the operation for 
each value of the grouping var stu_id, the values of the grouping var in BY 
statement are not transposed */
   var score1 - score3;
   /*values of score1-3 will be transposed*/
run;

/*This is almost what you want. All that is needed is to 
rename the variable COL1 to all_score, 
rename the _NAME_ , and eliminate _label_ variables, 
and remove the observations with missing score values. */

proc transpose data=oneper
               out=manyper_final (rename = (col1 = all_score _name_ = score_type)
               						drop =  _label_
               						where=(all_score ne .));
   by stu_id;
   var score1 - score3;
run;

proc print data = oneper;
title 'one obs per student';
run;
proc print data = manyper_final;
title 'multiple obs per student';
run;

proc sort data = oneper ;
by stu_id name gender;
run;
proc transpose data=oneper
               out=manyper_final_extra (rename = (col1 = all_score _name_ = score_type)
               						drop =  _label_
               						where=(all_score ne .));
   by stu_id name gender;
   var score1 - score3;
run;
proc print data = manyper_final_extra;
title 'multiple obs per student_extra vars';
run;
