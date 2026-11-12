/**create permanent SAS datasets using LIBNAME Statement**/

LIBNAME score "/folders/myfolders";
/* it refers to my local folder C:\SASUniversityEdition\myfolders 
(SAS studio is a visualization software, 
if you are like me using SAS studio, 
you need to create a shared folder ‘SASUniversityEdition\myfolders’ 
on your local drive so SAS studio can see and use your local files.)

If you have SAS installed on your local computer, the directory will be 
like this: 
LIBNAME  score "C:\your-folder"; */

DATA score.stu_math;
INPUT stu_name $ math_score;
Datalines;  		
Speedy 90
Tim 91
Sharon 100
Mike 60
;      
run;      

/* create a temparary SAS dataset   --- check the WORK library*/

DATA stu_math_temp;
INPUT stu_name $ math_score;
Datalines;  		
Speedy 90
Tim 91
Sharon 100
Mike 60
;      
run;        

/* close current SAS session, open a new SAS session, 
look at library WORK, Score*/
 
 LIBNAME  score "/folders/myfolders";
 
 proc print data = score.stu_math;
 run;
