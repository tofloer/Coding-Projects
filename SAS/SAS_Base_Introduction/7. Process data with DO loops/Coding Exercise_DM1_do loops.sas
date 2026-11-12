/*Section: Generating data with DO loops
coding exercise:

instruction:
1. compute each year's salary you will have 
if you start with $60,000 at a 3% increase rate per year for 5 years.
2. Conditionally Executing DO Loops:
how many years it will take the salary to reach $100,000 per year 
if the salary increase 3% each year, again the initial salary is $60,000.
*/

/*1*/
data salaryincrease1 (drop= counter);
   Increase_rate = .03;
   salary = 60000;
   do counter = 1 to 5;
      salary+salary*increase_rate;
      Year+1;
      output;
   end;
   format salary dollar10.2;
run;

/*2*/
data salaryincrease2;
salary = 60000;
   do until(salary>=100000);
      salary+salary*.03;
      Year+1;
   end;
run;