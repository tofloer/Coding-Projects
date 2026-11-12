/*proc sql: retrieve data from a single table using Select statement*/
proc import datafile = "/folders/myfolders/score_data_id" 
DBMS = xlsx out = score_data ;
run;

proc sql;
   create table scoredata0 as
   select stu_id,
          gender,
          name
   from score_data
   where gender in  ('m');
quit;

/*if you don't include the Create keyword 
you are not creating an output data set, so, 
by default, the result of this query is sent 
as a listing in the Output window */

proc sql;
   select stu_id,
          gender,
          name
   from score_data
   where gender in  ('m');
quit;

/*If you want to select all the variables from a data set, 
you can use an asterisk (*)*/
proc sql;
   select *
   from score_data
   where gender in  ('m');
quit;
