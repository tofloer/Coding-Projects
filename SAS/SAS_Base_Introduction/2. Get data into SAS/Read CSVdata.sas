/* Read .csv data*/

/* Can also look at the original datafile DATA_commas.csv 
5 obs and 4 vars

Tim, M, 50, 145
Sara, , 23, 130 
Mike, M, 65, 180
Laura, F, ,130
Sean, M, 15, 167
*/



/* For detailed explanation on INFILE and INPUT statements,
Please check my previous tutorial on 
how to read data seperated by blanks */
 
data sdata_commas;
   infile "/folders/myfolders/DATA_commas.csv" dsd;
   input name $ Gender $ Age Weight;
run;

proc print data=sdata_commas;
run;
