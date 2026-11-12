/* Enhance frequency tables with options (NLEVELS, ORDER = )*/

/*using data step program with If Then Else statements
to create data set to be used in Proc Freq to 
generate frequency tables --- var Gender and Grade will be used*/
proc import datafile = "/home/c.sharonx/score_data_miss" 
DBMS = xlsx out = scoredata0 replace ;
run;
data scoredata1;
   set scoredata0; 
   /*using SAS functions*/
    TotalScore = sum (score1, score2, score3);
	AverageScore = mean (score1, score2, score3);
run;

DATA scoredata_IfThenElse;
set scoredata1;
	/*using If...Then...Else statements*/
		/*gender_num*/
	If gender = 'm' then gender_num = 1;
	else if gender = 'f' then gender_num = 0;
	else gender_num = . ;
	
	/*take*/
	IF score1 NE . AND score2 NE .  AND score3 NE .  THEN take = 'complete';
	else take = 'incomplete';
	
	/*grade & pass*/
	IF AverageScore >= 90 THEN DO;
       grade = 'A'; 
	   pass = 'pass';
    end;
    Else If averageScore >= 80 then do;
       grade = 'B';
       pass = 'pass';
	END;  
	Else If averageScore >= 70 then do;
       grade = 'C';
       pass = 'pass';
	END;  
	Else If averageScore >= 60 then do;
       grade = 'D';
       pass = 'pass';
	END;  
	Else If 0 =< averageScore < 60 then do;
       grade = 'F';
       pass = 'fail';
	END; 
	else do;
	grade = ' ';
       pass = ' ';
	END; 
	   
run;    
 
/*options:
(1) NLEVELS option Displays the number of levels for all TABLES variables
By specifing the NLEVELS option in the PROC FREQ statement, 
PROC FREQ displays the "Number of Variable Levels" table.

(2) ORDER=	option
specifies the order of the variable levels in the frequency and crosstabulation tables, 
which you request in the TABLES statement.

ORDER=DATA | FORMATTED | FREQ | INTERNAL
The ORDER= option can take the following values:
DATA		Order of appearance in the input data set
FORMATTED	External formatted value, except for numeric variables with no explicit format, 
            which are sorted by their unformatted (internal) value
FREQ		Descending frequency count; levels with the most observations come first in the order
INTERNAL	Unformatted value

By default, ORDER=INTERNAL. The FORMATTED and INTERNAL orders are machine-dependent. 
The ORDER= option does not apply to missing values, which are always ordered first.*/
proc freq data = scoredata_ifthenelse nlevels  ; 
tables gender grade gender*grade/ nocol norow nopercent; 
title "frequency tables: nlevel";
run;

proc freq data = scoredata_ifthenelse nlevels order = FREQ; 
tables gender grade gender*grade /nocol norow nopercent; 
title "frequency tables: order = FREQ"; 
run; /*with order, highest freq category displayed first*/

/*(3). LIST option and CROSSLIST option in Table statement:
When three or more variables are specified, the multiple levels of n-way tables 
can produce considerable output. Such bulky, often complex crosstabulations are 
often easier to read when they are arranged as a continuous list. 

The LIST option
To generate list output for crosstabulations, add a slash (/) and 
the LIST option to the TABLES statement in your PROC FREQ step.
Note: this option eliminates row and column frequencies and percentages

The CROSSLIST option displays crosstabulation tables in ODS column format 
instead of the default crosstabulation cell format. 
In a CROSSLIST table display, the rows correspond to the crosstabulation table cells, 
and the columns correspond to descriptive statistics such as Frequency and Percent.
The CROSSLIST table displays the same information as the default crosstabulation table,
but uses an ODS column format instead of the table cell format */
proc freq data = scoredata_ifthenelse order = FREQ; 
tables gender*grade ; 
title "frequency tables: order = FREQ"; 
run;
proc freq data = scoredata_ifthenelse order = FREQ; 
tables gender*grade / list; 
title "frequency tables: list"; 
run;
proc freq data = scoredata_ifthenelse order = FREQ; 
tables gender*grade / crosslist; 
title "frequency tables: crosslist"; 
run;





