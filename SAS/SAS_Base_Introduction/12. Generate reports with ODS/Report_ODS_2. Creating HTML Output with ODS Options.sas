/*Creating HTML Output with ODS Options:

To create simple HTML output, you do not have to specify the ODS HTML statement. 
However, to create HTML output with options specified, you open the HTML destination 
using the ODS HTML statement.
________________________________________
General form, ODS HTML statement:
ODS HTML BODY=file-specification;
ODS HTML CLOSE;
----------------------------------------
There are diff ways for file-specification, We will use the complete directory path and HTML filename
in this SAS program.

In this eample, we will use the ODS HTML statement to direct the results from multiple procedures 
to the same HTML file.*/

proc import datafile = "/folders/myfolders/score_data_miss999" 
DBMS = xlsx out = scoredata0 replace ;
run;

ods html file = "/folders/myfolders/scoredata0.html" ;
proc print data=scoredata0;
id name;
title 'score data table';
run;
proc means data=scoredata0;
var score1 score2 score3;
class gender;
title 'score means by gender';
run;
ods html close;











