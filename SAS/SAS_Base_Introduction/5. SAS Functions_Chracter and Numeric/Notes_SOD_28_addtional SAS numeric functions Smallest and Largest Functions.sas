/*additional SAS numeric functions: SMALLEST and LARGEST

previously introduced (still included examples below):
sum : sum of a selection of variables
mean : mean of a selection of variables
max : the maximum values from a selection of variables
min: the minimum values from a selection of variables

SMALLEST(k, value-1 <, value-2...>)
Returns the kth smallest nonmissing value.
Required Arguments:
k: is a numeric constant, variable, or expression that specifies 
which value (i.e. kth smallest value) to return.
value: specifies a numeric constant, variable, or expression.

LARGEST(k, value-1 <, value-2 ...>)
Returns the kth largest nonmissing value.
Required Arguments:
k: is a numeric constant, variable, or expression that specifies 
which value (i.e. kth largest value) to return.
value: specifies the value of a numeric constant, variable, or expression to be processed.

If k is missing, less than zero, or greater than the number of values, 
the result is a missing value and _ERROR_ is set to 1. 
Otherwise, if k is greater than the number of nonmissing values, 
the result is a missing value, but _ERROR_ is not set to 1.
*/
proc import datafile = "/home/c.sharonx/score_data" 
DBMS = xlsx out = scoredata0 replace ;
run;

DATA scoredata3;
   set scoredata0; 
   
/*use SAS functions to use only non-missing values for the computation */    
	TotalScore_func = sum (score1, score2, score3);
	AverageScore_func = mean (score1, score2, score3);
	max_score = max (score1, score2, score3); 
	min_score = min (score1, score2, score3); 
	small_2nd_score = smallest (2,score1, score2, score3); /*returns 2nd smallest values of the three scores*/
	large_1st_score = largest (1,score1, score2, score3);/*returns 1st largest values of the three scores*/
RUN;
PROC PRINT DATA = scoredata3;
RUN;

