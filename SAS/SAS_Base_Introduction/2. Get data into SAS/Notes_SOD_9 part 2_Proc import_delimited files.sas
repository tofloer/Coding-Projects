/* Proc Import: Read delimited external data files:

DELIMITER Statement
Specifies the delimiter that separates columns of data in the input file.
Default:	
-Comma for .CSV files
-Blank space for all other file types

If you specify DBMS=DLM, you must also specify the DELIMITER= statement.
Note: If you omit DELIMITER=, the IMPORT procedure assumes that the delimiter is a space.

GUESSINGROWS Statement
Specifies the number of rows of the file to scan to determine the appropriate data type 
and length for the variables.
Default:	20
Keep in mind that the more rows you scan, the longer it takes for the PROC IMPORT to run

Syntax
GUESSINGROWS=n | MAX;
n: indicates the number of rows the IMPORT procedure scans in the input file to 
determine the appropriate data type and length of variables. The range is 1 to 2147483647 (or MAX). 
MAX:can be specified instead of 2147483647. Specifying the maximum value could adversely affect performance.

The GETNAMES= no for all codes since we don't have var-names in the first line of the data files
 */

/*the original blank-seperated input data

Tim M 50 145
Sara . 23 130 
Mike M 65 180
Laura F . 130
Sean M 15 167

*/
 
proc import datafile="/home/c.sharonx/DATA_blanks.txt"
   out=sdata_blanks
   dbms=dlm
   replace;
   GUESSINGROWS= 2; 
   *delimiter=' '; /*this may be omitted as it is the default*/
   getnames=no;
run;
proc print data=sdata_blanks;
run; /*extra variable with no data values*/
data sdata_blanks_1;
set sdata_blanks;
drop var5;
run;
proc print data=sdata_blanks_1;
run;


/* .csv file
Tim, M, 50, 145
Sara, , 23, 130 
Mike, M, 65, 180
Laura, F, ,130
Sean, M, 15, 167
*/


proc import datafile="/home/c.sharonx/DATA_commas.csv"
   out=sdata_comma
   dbms=dlm
   replace;
   delimiter=',';
   getnames=no;
run;
proc print data=sdata_comma;
run;

/*Again same data as used for blanks seperated or CSV data files,
only change the delimiter to colons :

Tim: M: 50: 145
Sara: : 23: 130 
Mike: M: 65: 180
Laura: F:  :130
Sean: M: 15: 167
*/
proc import datafile="/home/c.sharonx/other_del_data.txt"
   out=sdata_colon
   dbms=dlm
   replace;
   delimiter=':';
   getnames=no;
run;
proc print data=sdata_colon;
run; /*extra obs without data values*/
data sdata_colon_1;
set sdata_colon;
if var1 = ' ' then delete;
run;
proc print data=sdata_colon_1;
run;

