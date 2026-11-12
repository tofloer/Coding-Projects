/*constructing DO loops_part2*/

/*The purpose of this program: 
compute the total amount of money you will have 
if you start with $100 and invest it 
at a 3.75% interest rate for 3 years. */


/*Drop counter variable*/

data earnings3 (drop= counter);
   Interest = .0375;
   Total = 100;
   do COUNTER = 1 to 3;
      Total + Interest*Total;
      Year+1;
   end;
   format Total dollar10.2;
run;

title "Listing of earnings with DO loops";
proc print data=earnings3 ;
run;


/*program with DO loops_OUTPUT*/

data earnings2;
   Interest = .0375;
   Total = 100;
   do COUNTER = 1 to 3;
      Total + Interest*Total;
      Year+1;
      output;
   end;
   format Total dollar10.2;
run;

title "Listing of earnings with DO loops_every year";
proc print data=earnings2 ;
run;

/*program without DO loops_OUTPUT*/
data earnings1;
   Interest = .0375;
   Total = 100;

   Year + 1;
   Total + Interest*Total;
   output;

   Year + 1;
   Total + Interest*Total;
   output;

   Year + 1;
   Total + Interest*Total;
   output;

   format Total dollar10.2;
run;

title "Listing of earnings without DO loops_every year";
proc print data=earnings1;
run;

