/*section: creating labels and format

Coding Exercise 2:

data: 
sale.xlsx

Instruction:
1. Re-run the codes in coding Exercise 1
2. create format 'salegroup' for averagesale with following categories
and store the format in library 'myfmts'
low-<700 = 'need improvement'
700-<900 = 'good'
900-high = 'top sale'
3. in Proc Print: reference/use the format 'salegroup' stored in the Myfmts library
as if in a new SAS session, meaning need to include LIBNAME statement
*/

/*1. Re-run the codes in coding Exercise 1*/
proc import datafile = "/folders/myfolders/sale" 
DBMS = xlsx out = sale0 replace ;
run;

data sale1;
set sale0;
   label emid	= 'employee id'
		sale_m1 = 'sale in Jan.'
		sale_m2 = 'sale in Feb.'
		sale_m3 = 'sale in Mar.'
		averagesale = 'average sale of quarter1';
	averagesale = mean (of sale_m1, sale_m2, sale_m3);
	format averagesale sale_m1 --sale_m3 dollar10.2;
run; 

proc print data = sale1 label;
run;

/*2*/
libname myfmts "/folders/myfolders/formats" ;
PROC FORMAT library = myfmts;                            
   VALUE salegroup low-<700 = 'need improvement'
   					700-<900 = 'good'
   					900-high = 'top sale'  ;
run;

/*3*/
libname myfmts "/folders/myfolders/formats" ;
options fmtsearch=(myfmts work library);
proc print data = sale1;
	format averagesale salegroup.;
run;


