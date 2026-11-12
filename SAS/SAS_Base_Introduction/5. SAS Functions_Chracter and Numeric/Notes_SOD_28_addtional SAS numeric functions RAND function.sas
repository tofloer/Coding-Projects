/*additional SAS numeric functions: 
RAND function 
Generates pseudo-random numbers from a distribution that you specify.
RAND(' dist', parm, ...parm)
'dist': is a character constant, variable, or expression that identifies the distribution
		for example:  Chi-Square, F, Normal, Uniform, etc.
There are many distributions, for complete information: 
https://documentation.sas.com/doc/en/vdmmlcdc/8.1/ds2ref
/p11xdsmtnx21pnn11ewv8nnblv4u.htm#p0puykspkij2ssn1ml76oddwodej
parm: are shape, location, or scale parameters that are appropriate for the specific distribution.

/*1. To generate random numbers between a and b, we can do so easily in the following SAS code
x= RAND('UNIFORM', <a, b>)
Range:	
The range of x is a ≤ x ≤b. 
The parameters a and b: default values are a=0 and b=1 -- without specifying a, b 
You must specify values for both a and b if you do not want to use the default values.*/

data k1;
	do i = 1 to 10;
		rand_num1= rand("Uniform"); /*default a=0, b=1*/
		rand_num2= rand("Uniform", 1,10);
		output;
	end;
run; 

/*2. Generating Random Integers:
If you have access to newer version of SAS, you can generate random integers directly.
The RAND function supports the "Integer" distribution and supports parameters 
that specify the min and max. 
To generate a random integer between 1 and 10 you can use the following statements:*/
data k2;
    do i = 1 to 10;
        rand_int = rand("integer",1,10);
        output;
    end;
run;

/*Reproducing a Random Number Stream
If you want to create reproducible streams of random numbers, 
use the CALL STREAMINIT routine to specify a seed value for random number generation, 
we need to use the CALL STREAMINIT routine before any invocation of the RAND function. */

/*The following code shows how to re-generate a variable in SAS that contains 
10 random values between 1 and 20
without CALL STREAMINIT routine, RAND will generate diff randon values at each run */
data k3;
   call streaminit(7); /*streaminit() function to ensure that this example is reproducible. */
   do i = 1 to 10;
   x = rand("integer", 1, 20);
   output;
   end;
run;




