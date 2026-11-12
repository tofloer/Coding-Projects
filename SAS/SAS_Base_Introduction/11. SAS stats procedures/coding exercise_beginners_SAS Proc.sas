/* SAS programming for beginners:
section: SAS statistical procedures

coding exercise
1. Use data set ‘scoredata1’ from the coding exercise for ‘data preparation’ (for details Please see step 1 and 2 in the coding exercise for ‘data preparation’)

2. Sort data by Gender, Averagecore; Print the sorted data set

3. Generate frequency tables for all character variables excluding Name

4. Generate statistic outputs for Ns1,2,3 and Averagecore by Gender using Proc Means and Proc Univariate
*/

proc import datafile = "/folders/myfolders/score_data_miss777" 
DBMS = xlsx out = scoredata0 replace ;
run;

data scoredata1;
set scoredata0;

ARRAY sc (3) score1 score2 score3;    
   ARRAY new (3) ns1 ns2 ns3; 
   DO i = 1 TO 3;                       
      IF sc(i) = 777 THEN new(i) =.;   
      Else if sc(i) NE 777 then new(i) = sc(i);
   END;  
   
averagescore = mean (ns1, ns2, ns3);
run;  

PROC SORT DATA = scoredata1 OUT = scoredata2 ;
by Gender descending averagescore;
run;

proc print data = scoredata2;
run;

proc freq data = scoredata2;
table gender;
run;

proc means data = scoredata2;
vars ns1 - ns3;
by gender;
run;

proc univariate data = scoredata2;
var ns1 - ns3;
by gender;
run;


