data stu_math;
input stu_name $ math_score;
cards;
Speedy 90
Tim 91
Sharon 100
Mike 60
Michael 75 
;
run;

proc means data = stu_math;
var math_score;
run; 

