/*Error Handling: identify & correct errors: 
The PUTLOG Statement
The Basics of Logic Errors:
A logic error occurs when the program statements execute, 
but produce incorrect results. Because no notes are written to the log, 
logic errors are often difficult to detect. 
Use the PUTLOG statement in the DATA step to help identify logic errors.

In PULOG, we will use the named output to display/write the value of variables, or formatted values, 
or to write values of all variables

the named output: list the variable name followed by an equal sign.
For example, this PUTLOG statement

putlog name= age=;
writes the values for NAME and AGE to the SAS log: 

name=Peterson age=21
name=Morgan age=17
*/
proc import datafile = "/folders/myfolders/score_data" 
DBMS = xlsx out = scoredata0 replace ;
run; /*look at data*/

/*Example 1: The DATA Step Produces Wrong Results 
but There Are No Error Messages*/
DATA scoredata1;
   set scoredata0;                                           
   AverageScore = mean (score1 + score2 + score3);      
PUTLOG 'the score infomation: ' Name= Score1= Score2= Score3= AverageScore= 5.2;
/*Use PUTLOG to write the value of a variable, formatted values, or to write values of all
variables*/
   If averagescore <60; /*with Conditional logic*/
RUN;

/*correct and then resubmit*/
DATA scoredata1_c;
   set scoredata0;                                           
   AverageScore = mean (score1 , score2 , score3);      
   If averagescore <60; /*with Conditional logic*/
   PUTLOG 'WARNING: Score below the minimum: ' Name= Score1= Score2= Score3= AverageScore= 5.2; 
/*PUTLOG: with Conditional logic*/
RUN;

/* Use temporary variables N and ERROR to debug a DATA step:
_N_ 
The number of times the DATA step iterated
Displays debugging messages for a specified number of iterations of the DATA step
_ERROR_
Initialized to 0, set to 1 when an error occurs
Displays debugging messages when an error occurs 
ERROR statement
The ERROR statement sets the automatic variable _ERROR_ to 1. 
Writing a message that you specify to the SAS log is optional. 
When _ERROR_ = 1, SAS writes the data lines that correspond to 
the current observation in the SAS log.*/
DATA scoredata1_e;
   set scoredata0; 
   AverageScore = mean (score1 + score2 + score3);   
   Error '???' Name= Score1= Score2= Score3= AverageScore= 5.2; 
   If averagescore <60; /*with Conditional logic*/
RUN;






