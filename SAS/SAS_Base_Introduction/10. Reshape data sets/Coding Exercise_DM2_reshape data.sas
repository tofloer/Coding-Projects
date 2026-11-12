/*section: Reconstruct/Reshape SAS Data sets in DATA step 
and using Proc TRANSPOSE

Coding Exercise:
data:
Weight_loss.xlsx
instruction:
0. change value 9999 to missing (periods)
1. transform data from one record per PID (patient ID)
to multiple record per PID using Proc Transpose
2. Values of 3 weight vars (weight0 weight1 weight2)
will be transposed 
3. all vars in the input data will be included in the transposed data
*/

proc import datafile = "/folders/myfolders/Weight_loss" 
DBMS = xlsx out = wl0 replace ;
run;
data wl1;
set wl0;
array v(2)	weight1 weight2;
   DO i = 1 TO 2;                       
      IF v(i) = 9999  THEN v(i) =.;   
   END; 
drop i;
run;

proc sort data = wl1 ;
by pid gender walk_steps;
run;
proc transpose data=Wl1
               out=WL_m (rename = (col1 = all_weight _name_ = weight_time)
               						drop =  _label_
               						where=(all_weight ne .));
   by pid gender walk_steps;
   var weight0-weight2;
run;
proc print data = wl_m;
title 'weight loss: multiple records per patient id';
run;
