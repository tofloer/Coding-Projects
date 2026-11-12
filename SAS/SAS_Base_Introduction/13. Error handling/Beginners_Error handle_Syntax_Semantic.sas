/*Error Handling: identify & correct errors: 
Types of Errors
SAS can detect several types of errors. Here are two common ones:
•	syntax errors that occur when program statements do not conform 
to the rules of the SAS language, such as spelling mistake, missing RUN statement,
missing semicolon,  or unbalanced quotation mark

•	semantic errors that occur when you specify a language element that is 
not valid for a particular usage (invalid option)

When a program that contains an error is submitted, messages about the error appear 
in the Log window. Here is what SAS does:
•	displays the word ERROR
•	identifies the possible location of the error
•	gives an explanation of the error

In SAS studio:
1.To correct errors in programs, edit them in the code editing window. 
2. Resubmitting a Revised Program 
(SAS studio clears the messages from the previous submission in the Log window )
(If you are working in the SAS windowing environment, it is a good idea 
to clear the messages from the Log window before resubmitting the program 
so that you do not confuse the old error messages with the new messages.)
3. Remember to check the Log window again to verify that your program ran correctly.

*/
proc import datafile = "/folders/myfolders/score_data" 
DBMS = xlsx out = scoredata0 replace  ;
run;

DATA scoredata1_c;
   set scoredata0;                                           
   AverageScore = mean (score1 , score2 , score3);      
   If averagescore <60;
run;

/*syntax errors*/
/*misspelled keyword*/
Proc mean data = scoredata1_c;
var averagescore;
run;

/*unmatched quotation marks

Sometimes a warning appears in the SAS log that indicates the following:
•	A quoted string has become too long.
•	A statement that contains quotation marks (such as a TITLE or FOOTNOTE statement) 
is ambiguous because of invalid options or unquoted text.
Sometimes, you don't even see warnings in LOG 

To correct, in SAS Studio, simply add the balancing quotation mark 
and resubmit the program. 
(However, in some environments, you may need to resolve the unbalanced quotation mark 
by canceling the submitted statements)*/
Proc means data = scoredata1_c;
var averagescore;
title 'info on average score'; 
/*need to delete the quotation mark after score to create a unbalanced quotation mark error*/
run;

/*missing semicolons ; 
When you omit a semicolon, SAS reads the statement that lacks the semicolon 
(along with the following statement) as one long statement.

-- Find the statement that lacks a semicolon. You can usually find it 
by looking at the underscored keywords in the error message and working backward*/
DATA scoredata1_c;
   set scoredata0;                                           
   AverageScore = mean (score1 , score2 , score3);      
   If averagescore <60
run;

/*semantic errors:
Invalid option
An invalid option error occurs when you specify an option 
that is not valid in a particular statement. 

When a SAS statement that contains an invalid option is submitted, a message 
appears in the Log window indicating that the option is not valid or not recognized.
*/

Proc means data = scoredata1_c label;
var averagescore;
run;
