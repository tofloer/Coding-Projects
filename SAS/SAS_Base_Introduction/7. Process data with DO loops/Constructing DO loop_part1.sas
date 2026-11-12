/*constructing DO loops_part1*/

/*The purpose of this program: 
compute the total amount of money you will have 
if you start with $100 and invest it 
at a 3.75% interest rate for 3 years. */

/*program without DO loops*/
data earnings1;
   Interest = .0375;
   Total = 100;

   Year + 1;
   Total + Interest*Total;
   

   Year + 1;
   Total + Interest*Total;
  

   Year + 1;
   Total + Interest*Total;
  

   format Total dollar10.2;
run;

title "Listing of earnings without DO loops";
proc print data=earnings1;
run;

/*program with DO loops*/

data earnings2;
   Interest = .0375;
   Total = 100;
   do Year = 1 to 3;
      Total + Interest*Total;
   end;
   format Total dollar10.2;
run;

title "Listing of earnings with DO loops";
proc print data=earnings2 ;
run;

/*After the third execution of the DO loop, 
the value of year is incremented to 4. 
Because 4 exceeds the stop value of the iterative DO statement, 
the DO loop stops executing, and processing continues to the 
next DATA step statement. The end of the DATA step is reached, 
the values are written to the Earnings data set, 
year shows as 4 instead of 3
and only one observation is written to the data set.
*/