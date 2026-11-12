/*Section: More on Arrays
coding exercise 1:

data:
allscore_miss_text.xlsx

instruction:
using Array to replacing 'missing' with blanks in vars
read math science write
then convert the character vars read math science write
to numeric vars
*/

proc import datafile = "/folders/myfolders/allscore_miss_text" 
DBMS = xlsx out = as0 replace ;
run;

data as1 (drop = i);
set as0;
array oldvar (4) read--write;
array newvar (4) readN mathN scienceN writeN;
do i = 1 to 4;
if oldvar(i) = 'missing' then oldvar(i) = ' ';
newvar(i) = input (oldvar(i), 8.);
end;
run;

  

