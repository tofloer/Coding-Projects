/*Conditionally Executing DO Loops*/

/*want to know how many years it will take to earn $50,000 
if you deposit $2,000 each year into an account that earns 10% interest. 
The DATA step below uses a DO UNTIL statement to perform 
the calculation until $50,000 is reached. Each iteration of 
the DO loop represents one year.*/


/*DO UNTIL*/
data invest;
   do until(Capital>=50000);
      capital+2000;
      capital+capital*.10;
      Year+1;
   end;
run;
/*During each iteration of the DO loop,
	2000 is added to the value of Capital to reflect the annual deposit 
of $2,000
	10% interest is added to Capital
	the value of Year is incremented by 1.
Because there is no index variable in the DO UNTIL statement, 
the variable Year is created in a sum statement 
to count the number of iterations of the DO loop. 
This program produces a data set that contains the single observation. 
To accumulate more than $50,000 in capital requires 13 years 
(and 13 iterations of the DO loop). */

/*DO WHILE*/
data invest;
   do while(Capital>=50000);
      capital+2000;
      capital+capital*.10;
      Year+1;
   end;
run;
/*in this program, because the value of Capital is initially zero, 
which is less than50,000, the DO loop does not execute.*/
