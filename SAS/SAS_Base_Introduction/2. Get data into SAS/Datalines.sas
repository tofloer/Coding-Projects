/*  Datalines   */

DATA stu_math;
INPUT stu_name $ math_score;
Datalines;  		/*CARDS*/
Speedy 90
Tim 91
Sharon 100
Mike 60
;              /*The ; must be placed at a new line*/
RUN;

Proc print data = stu_math; 
Run;

