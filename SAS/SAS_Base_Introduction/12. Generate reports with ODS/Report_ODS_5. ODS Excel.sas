/* ODS Excel destination:
Customizing Your Excel Output:
create multiple Excel sheets according to by group variable's values*/

proc import datafile = "/folders/myfolders/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;

PROC SORT DATA = scoredata0 OUT = scoredata1 ;
by Gender ; 
run;

/*The SHEET_INTERVAL= option specifies that a new sheet is created after each by group. 
The SHEET_LABEL= option specifies the text “Gender =” as the first part of the worksheet label. 
The EMBEDDED_TITLES= option specifies that the title created by the TITLE statement is embedded 
in the Excel worksheet. 
The EMBED_TITLES_ONCE option specifies that the title is embedded once, at the top of each sheet.

The ODS NOPROCTITLE statement suppresses the writing of the procedure title that produces the results.*/

ods excel file="/folders/myfolders/multi_table.xlsx" 
   options(sheet_interval="bygroup"
           sheet_label="Gender ="
           embedded_titles="yes"
           embed_titles_once="yes" );
TITLE 'Summary of Scores by Gender';
ods noproctitle;
PROC MEANS data = scoredata1 MAXDEC = 1 n mean max min;      
   BY gender;                             
   VAR score1 score2 score3;      
   where gender ne ' ';
RUN;
ods excel close;

