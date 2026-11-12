/*PROC Contents: 
1. Viewing the Contents of SAS Libraries

You can use a PROC CONTENTS step to view the contents of a SAS library. 
When you specify the keyword _ALL_ in the PROC CONTENTS statement, 
the step displays a list of all the SAS files that are in the specified 
SAS library.
By default, a PROC CONTENTS report includes the descriptor portion of each data set 
in the SAS library (The descriptor portion of a SAS data set contains information 
about the attributes of the data set and information about each
variable.) . You can use the NODS option to suppress the descriptor portions 
in the report*/
LIBNAME score "/folders/myfolders"; 
/*After you end your SAS session, the contents of a permanent library still exist 
in their physical location. However, each time you start a new SAS session, 
you must resubmit a LIBNAME statement to reassign the libref to that permanent library*/
PROC CONTENTS DATA=score._ALL_ NODS; 
RUN;
PROC CONTENTS DATA=score._ALL_ ; 
RUN; 

/*2. Viewing the contents of SAS Data Sets*/


PROC CONTENTS DATA=score.score_data;
RUN;
