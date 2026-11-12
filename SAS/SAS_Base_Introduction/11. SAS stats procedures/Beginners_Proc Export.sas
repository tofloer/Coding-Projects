/*Proc EXPORT: 
Export data to create comma-delimited and Excel data files.

Procedure features:	
The EXPORT procedure statement arguments:
DATA=
DBMS=
OUTFILE=
REPLACE
*/

/* Example: Exporting a Subset of Observations to a CSV File
The WHERE option requests a subset of the observations.
The OUTFILE option specifies the output file. 
The DBMS option specifies that the output file is a CSV file, 
The REPLACE option overwrites the target CSV, if it exists.*/

LIBNAME score "/folders/myfolders";
/*After you end your SAS session, the contents of a permanent library still exist 
in their physical location. However, each time you start a new SAS session, 
you must resubmit a LIBNAME statement to reassign the libref to that permanent library*/
proc export data=score.score_data (where=(gender='f'))
   outfile='/folders/myfolders/Femalelist.csv'
   dbms=csv
   replace;
run;/*then view the file as a text by right click the data file*/

/* The DBMS= option is dlm, which
means the output file is a delimited file. */
proc export data=score.score_data
outfile="/folders/myfolders/sd"
dbms=dlm
replace;
delimiter='&';
run;  

/* Example: Exporting to a Excel file
The data are written into a worksheet called 'data'. */
proc export data=score.score_data
outfile='/folders/myfolders/sd_excel'
dbms=XLSX replace;
sheet='data';
run; 
