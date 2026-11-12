/*FIND function:
The FIND function enables you to search for a specific substring of characters 
within a specified character string.
--- If the substring is not found in the string, the FIND function returns a value of 0. 
--- If there are multiple occurrences of the substring, the FIND function returns only 
the position of the first occurrence.

syntax: FIND(string,substring<,modifiers><,startpos> ) 

The FIND function and the INDEX function both search for substrings of characters 
in a character string. However, the INDEX function does not have the modifier 
nor the start-position arguments.*/

proc import datafile = "/home/c.sharonx/Chara_data1" 
DBMS = xlsx out = scoredata0 replace ;
run;

/*FIND(string,substring<,modifiers><,startpos> ) 
required arguments:
(1) string specifies a character constant, variable, or expression that is searched 
for substrings. 
(2) substring is a character constant, variable, or expression that specifies the 
substring of characters to search for in string. 
Note: If string or substring is a character literal, you must enclose it in quotation marks.*/
data scoredata1;
set scoredata0;
Find_m = find(name,'m');
run;
proc print data = scoredata1;
var name find_m;
title "find(name,'m')";
run;
/*optional arguments:
(3) modifiers is a character constant, variable, or expression that specifies 
one or more modifiers. 
---The modifier i causes the FIND function to ignore character case during the search. 
If this modifier is not specified, FIND searches for character substrings with 
the same case as the characters in substring.
---The modifier t trims trailing blanks from string and substring.
--- Note: If modifier is a constant, enclose it in quotation marks.
Specify multiple constants in a single set of quotation marks. 
Modifier can also be expressed as a variable or an expression.*/
data scoredata2;
set scoredata0;
Find_m = find(name,'m','i');
run;
proc print data = scoredata2;
var name find_m;
title "find(name,'m', 'i')";
run;
/*(4) startpos is an integer that specifies the position at which the 
search should start and the direction of the search. 
The default value for startpos is 1. 
--- If start-position is not specified, FIND starts the search at the beginning of 
the string and searches the string from left to right. 
--- If start-position is specified, the absolute value of start-position determines 
the position at which to start the search. 
--- The sign of start-position determines the direction of the search.
when startpos is positive, starts the search at position start-position,
FIND searches from startpos to the right,
When startpos is negative, starts the search at position start-position,
FIND searches from startpos to the left.*/
data scoredata3;
set scoredata0;
Find_m = find(name,'m','i',2);
run;
proc print data = scoredata3;
var name find_m;
title "find(name,'m','i',2)";
run;
data scoredata4;
set scoredata0;
Find_m = find(name,'m','i',-1);
run; /*start position is 1, search at position 1 to the left*/
proc print data = scoredata4;
var name find_m;
title "find(name,'m','i',-1)";
run;
 



