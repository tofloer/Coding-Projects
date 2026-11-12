/* describing data using Proc Univariate

PROC UNIVARIATE is a powerful Procedure that combines many of the features found
in other analytical PROC such as FREQ, MEANS,
SUMMARY, and TABULATE into a single PROC step.
Proc UNIVARIATE is an excellent exploratory
data analysis tool. In some cases it provides information that cannot be
found on any other BASE SAS PROC, such as
information on the data's median, mode, quartiles
and percentiles. */

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;  

PROC UNIVARIATE data = scoredata1;                                   
   VAR averagescore;      
   TITLE 'Summary of Average scores';
RUN;

/*The BY Statement
The BY Statement on PROC UNIVARIATE allows
us to obtain separate sub-group analyses for each
value of the BY variable.
Whenever using the BY statement it requires that
the data be in the BY order. */

PROC SORT DATA = scoredata1 OUT = scoredata2 ;
by Gender ;
run;
PROC UNIVARIATE data = scoredata2;                                   
   VAR averagescore;   
   BY gender;
   TITLE 'Summary of Average scores by gender';
RUN;




