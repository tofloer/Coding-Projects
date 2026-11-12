/*inner join::: 
It is equivalent to a merge in a DATA step where only the obs from both data sets 
Are kept. you don’t have to sort the data sets first when you use SQL*/

proc import datafile = "/folders/myfolders/score_data_id_no_gender_1l" 
DBMS = xlsx out = score_data_ng ;
run; /*stu_id = 5 in not included*/

proc import datafile = "/folders/myfolders/score_data_id_gender_only_1l" 
DBMS = xlsx out = score_data_g ;
run; /*stu_id = 11 in not included*/

/*Using keyword Join*/
proc sql;
   CREATE TABLE inner_sql AS
   SELECT g.stu_id, 
          
          score1, 
          score2, 
          score3,
          gender, g.Name
/*For same-named variables in two tables, prefix the variable name with 
the table name or aliases*/       
   FROM score_data_g AS g INNER JOIN  score_data_ng AS ng
/*To make it easier to name variables with the same name from different 
tables, you create aliases for each of the tables*/
   ON score_data_g.stu_id =  score_data_ng.stu_id;
/*when you use the keyword JOIN to join two 
tables, you use an ON clause instead of a WHERE clause. */
        
quit;

/*using Where clause*/
proc sql;
   CREATE TABLE inner_sql_where AS
   SELECT g.stu_id, 
          
          score1, 
          score2, 
          score3,
          gender, g.Name     
   FROM score_data_g AS g,  score_data_ng AS ng
   where score_data_g.stu_id =  score_data_ng.stu_id;
quit;

proc print data=score_data_g;
title "Students gender information - stu_id = 11 in not included";
run;
proc print data=score_data_ng;
title "Students scores information - stu_id = 5 in not included";
run;

proc print data=inner_sql;
title "Performing an Inner Join Using Proc SQL using Join";
run;

proc print data=inner_sql_where;
title "Performing an Inner Join Using Proc SQL using Where clause";
run;


/*merge in a Data step*/
proc sort data = score_data_g;
by stu_id;
run;

proc sort data = score_data_ng;
by stu_id;
run;

data MMA;
merge score_data_g(in = A) score_data_ng (in = B);
/*Creating Temporary IN= Variables a and b*/
/*if score_data_g contributed data to the current observation,
  A = 1, otherwise A = 0, etc*/
by stu_id;
if A = 1 and B = 1;
run;

proc print data=mma;
title "a merge in a data step";
run;
proc print data=inner_sql;
title "Performing an Inner Join Using Proc SQL using Join";
run;





