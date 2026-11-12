/* SAS Macro Essentials
coding exercise:
data: 
allscore.xlsx
instruction: 
generate proc univariate outputs for all 4 score variables 
1. using macro variable replace score vars' names
2. using macro with parameter
 */

proc import datafile = "/folders/myfolders/allscore" 
DBMS = xlsx out = score0 replace ;
run;

/*1. using macro variable replace score vars' names*/
%let score = read;  /*run %let code for each score var one at a time*/
*%let score = math;  
*%let score = science;
*%let score = write;

/*The macro processor resolves the reference to the macro variable SCORE*/
proc univariate data = score0 ;
var &score;
title "Proc Univariate output for &score";
run;

/*2. using macro with parameter*/

%macro score_info(score_var= );
   proc univariate data = score0;
      var &score_var;
   run;
%mend score_info;

%score_info(score_var= read)
%score_info(score_var= write)
%score_info(score_var= math)
%score_info(score_var= science)
