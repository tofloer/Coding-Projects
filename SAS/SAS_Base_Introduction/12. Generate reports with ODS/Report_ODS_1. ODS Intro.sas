/* ODS Introduction 

The SAS Output Delivery System (ODS) gives you greater flexibility in generating, storing, and reproducing 
SAS procedure and DATA step output along with a wide range of formatting options. ODS provides formatting 
functionality that is not available when using individual procedures or the DATA step without ODS enhancements.

You can produce different types of formatted output with ODS:::
EPUB Output: With the ODS EPUB statement, you can create an e-book for use with e-book readers.
Excel Output: With ODS, you can produce tabular output, which can be viewed with Excel.
HTML Output: With ODS, you can produce output in HTML 
PDF Output: With ODS, you can produce output in PDF 
RTF Output: With ODS, you can produce RTF (Rich Text Format) output, which is used with Microsoft Word.
PowerPoint Output, etc.

An ODS block has the following form:
________________________________________

ODS output-destination 1 <options(s)>;
...
ODS output-destination (n) <options(s)>
<your SAS program>
ODS destination close statement 1;
...
ODS destination close statement (n)
________________________________________
In the ODS block, output-destination is the name of a valid ODS destination and option(s) are options that 
are valid for that destination. Your SAS program is inserted between the beginning ODS destination statement 
and the ODS CLOSE statement.
The FILE= or BODY= option is valid for most ODS statements. FILE= is an alias for BODY=.
These options specify the name and path of the file that you are generating.

HTML output is the default destination. Therefore, when you request another destination, your programs create 
both HTML output and output in the new destination. To prevent HTML output from being created, use this statement:

ods html close;
*/

proc import datafile = "/folders/myfolders/score_data_miss999" 
DBMS = xlsx out = scoredata0 replace ;
run;

ods pdf file="/folders/myfolders/scoredata0.pdf" ;
ods rtf file="/folders/myfolders/scoredata0.rtf" ;
ods powerpoint file = "/folders/myfolders/scoredata0.ppt" ;
proc print data=scoredata0;
id name;
run;
ods pdf close;
ods rtf close;
ods powerpoint close;
*ods _ALL_ CLOSE;

/***open outputs from the file diectory*/


