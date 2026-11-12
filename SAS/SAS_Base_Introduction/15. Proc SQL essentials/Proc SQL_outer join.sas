/* outer joins:::
left, right and full join*/

proc import datafile = "/folders/myfolders/score_data_id_no_gender_1l" 
DBMS = xlsx out = score_data_ng ;
run; /*no stu_id = 5*/

proc import datafile = "/folders/myfolders/score_data_id_gender_only_1l" 
DBMS = xlsx out = score_data_g ;
run; /*no stu_id = 11*/

proc print data=score_data_g;
title "Students gender information - stu_id = 11 in not included";
run;
proc print data=score_data_ng;
title "Students scores information - stu_id = 5 in not included";
run;
/*A left outer join lists matching rows and rows from the left-hand table 
(the first table listed in the FROM clause) that do not match any row 
in the right-hand table.*/ 
title 'left join';
proc sql;
   SELECT g.stu_id,ng.stu_id,
          
          score1, 
          score2, 
          score3,
          gender, g.Name
         
   FROM score_data_g AS g LEFT JOIN  score_data_ng AS ng
   ON score_data_g.stu_id =  score_data_ng.stu_id;
        
quit;

/*A right join: lists matching rows and rows from the right-hand table (the second 
table listed in the FROM clause) that do not match any row in the left-hand table. */
title 'right join';
proc sql;
   SELECT g.stu_id,ng.stu_id,
          
          score1, 
          score2, 
          score3,
          gender, g.Name
         
   FROM score_data_g AS g RIGHT JOIN  score_data_ng AS ng
   ON score_data_g.stu_id =  score_data_ng.stu_id;
        
quit;

/*A full join includes all rows from both tables*/
title 'full join';
proc sql;
   SELECT g.stu_id,ng.stu_id,
          
          score1, 
          score2, 
          score3,
          gender, g.Name
         
   FROM score_data_g AS g FULL JOIN  score_data_ng AS ng
   ON score_data_g.stu_id =  score_data_ng.stu_id;
        
quit;



