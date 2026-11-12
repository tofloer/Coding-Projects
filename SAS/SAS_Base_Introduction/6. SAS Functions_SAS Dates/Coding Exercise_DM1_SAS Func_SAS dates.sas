/*section: work with SAS dates 
Coding Exercise:
data:
patient_HD.xlsx
instruction:
1. extract the year of patients administered to hospital
2. calculate the number of days that patients stayed in hospital in two ways
3. creating a var containing today's date as analysis date 
*/

proc import datafile = "/folders/myfolders/patient_hd" 
DBMS = xlsx out = p0 replace ;
run;

data p1;
set p0;
h_year = year (start_date);

days_inH1=intck('day',start_date,end_date);
days_inH2 = datDIF(start_date,end_date,'ACT/ACT');


analysis_date = today();
format analysis_date date10.; /*you may use other format*/
run;

