/* Proc format CNTLIN= option 
CNTLIN=input-control-SAS-data-set
specifies a SAS data set from which PROC FORMAT builds informats or formats.
This program will show how to create a format from an input data set;
then use the created format in another data set 

There are three variables that are required in the formatting dataset. They are:
FMTNAME:  The format name
START:  The left side of the formatting = sign 
LABEL:  The right side of the formatting = sign

These three variables in a dataset will define a format that can also be defined using more traditional
code --- as we saw in previous tutorials:
PROC FORMAT;
 VALUE $fmtname
 ‘Start’ = ‘label’
 …etc….
 ; */

/*Create the variables Fmtname. The RETAIN statement is more efficient than an assignment statement in this case.
RETAIN retains the value of Fmtname in the program. Fmtname specifies the name Statefmt, which is the character format 
that the data set creates.*/
DATA fmt;
 Retain fmtname '$statefmt';
 Length start $2 label $2;
 Input start $ label $;
 datalines;
01 CA
02 TX
03 NJ
;
RUN;
proc print data = fmt; run;
PROC FORMAT CNTLIN=fmt;
RUN; 

proc import datafile = "/folders/myfolders/proc format_CNTLIN" 
DBMS = xlsx out = state0 replace ;
run;
data state1;
set state0;
format state $statefmt.;
run;



