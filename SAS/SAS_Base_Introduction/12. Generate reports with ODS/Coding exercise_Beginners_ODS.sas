/*ODS
Coding Exercise:
data:
sale_by_state.xlsx
instruction:
1. create Excel report using ODS
2. each Excel sheet is created for each state 
3. report sale average/mean for 3 sale variables sale 1- 3
*/
proc import datafile = "/folders/myfolders/sale_by_state" 
DBMS = xlsx out = ss replace ;
run;

proc sort data = ss; by state; run;
ods excel file="/home/c.sharonx/Sale Report by State.xlsx" 
   options(sheet_interval="bygroup"
           sheet_label="State ="
           embedded_titles="yes"
           embed_titles_once="yes" );
TITLE 'Sale Report by State';
ods noproctitle;
PROC MEANS data = ss MAXDEC = 1 n mean max min;      
   by state;                             
   VAR sale1 sale2 sale3;      
RUN;
ods excel close;