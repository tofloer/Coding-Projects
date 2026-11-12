*Read data seperated by blanks into SAS;

/* 
the original data

Tim M 50 145
Sara . 23 130 
Mike M 65 180
Laura F . 130
Sean M 15 167
*/



data sdata_blanks;
   infile "/folders/myfolders/DATA_blanks.txt";
 *C:\SASUniversityEdition\myfolders\DATA_blanks.txt’;
   input name $ Gender $ Age Weight;
run;

proc print data=sdata_blanks;
run;

