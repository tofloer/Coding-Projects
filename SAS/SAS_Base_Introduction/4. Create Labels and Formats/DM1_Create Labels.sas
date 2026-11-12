/*create label*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

/*label in data step*/
data scoredata1;
set scoredata0;
   label Score1 = 'Math Score'
         Score2 = 'Science Score'
         Score3 = 'English Score';
run;   

proc print data = scoredata1 label;
run;
/*use the LABEL option in the PROC PRINT statement to 
display labels*/

/*label in Proc step*/
proc print data = scoredata1 label;
label name = 'student name';  
/*The LABEL statement applies only to the PROC step in which it appears.*/
run;

/*Notice that the PROC PRINT statement still requires the LABEL option 
in order to display labels. 
Other SAS procedures display  assigned labels without additional option.*/
proc means data = scoredata1;
var score1;
run;

/*an example of using lable and split options in Proc print
LABEL
specifies to use the variables' labels as column headings.
SPLIT='split-character'
specifies the split character, which controls line breaks in column headings.
In this case, SPLIT= identifies the asterisk as the character that starts a new line in column headings.*/
proc print data = scoredata1 label split = '*' ;
label name = 'student * name';  
/*The LABEL statement applies only to the PROC step in which it appears.*/
run;












