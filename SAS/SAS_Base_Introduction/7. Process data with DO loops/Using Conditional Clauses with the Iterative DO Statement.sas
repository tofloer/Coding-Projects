/*Using Conditional Clauses with the Iterative DO Statement*/

/*In this DATA step, 
Suppose I want to limit the number of years I 
invest my capital to 10 years. And Add the UNTIL expression
to deternime years it takes for an investment to reach $50,000)
to further control the number of iterations. 

This iterative DO statement enables you to execute the DO loop 
until Capital is greater than or equal to 50000 or until 
the DO loop executes ten times, whichever occurs first.*/

data invest;
   do year=1 to 10 until(Capital>=50000);
      capital+2000;
      capital+capital*.10;
   end;
   if year=11 then year=10;
run;
/*In this case, the DO loop stops executing after ten iterations, 
and the value of Capital never reaches 50000. */


data invest;
   do year=1 to 10 until(Capital>=50000);
      capital+4000;
      capital+capital*.10;
   end;
   if year=11 then year=10;
run;

/*If you increase the amount added to Capital each year to 4000,
 the DO loop stops executing after the eighth iteration 
 when the value of Capital exceeds 50000.*/
