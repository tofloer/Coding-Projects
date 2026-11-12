/* Proc Import Excel file   */

/*no optional statements*/
proc import datafile = "/folders/myfolders/excel_data" 
DBMS = xlsx out = excel_data1 replace ;
run;

/*optional statements*/
/*getnames= No*/
proc import datafile = "/folders/myfolders/excel_data" 
DBMS = xlsx out = excel_data2 replace ;
getnames = no;  
/*By default, the IMPORT procedure will 
take the variable names from the first row of the spreadsheet.*/
run;

/* Sheet= "sheet-name" ; */
proc import datafile = "/folders/myfolders/excel_data" 
DBMS = xlsx out = excel_data3 replace ;
SHEET = "name_class";  /*bring in the first sheet by defaut if not specify*/
run;

/*range = "sheet-name$UL:LR" ; */
proc import datafile = "/folders/myfolders/excel_data" 
DBMS = xlsx out = excel_data4 replace ;
RANGE = "name_class$A1:B6"; 
run;








