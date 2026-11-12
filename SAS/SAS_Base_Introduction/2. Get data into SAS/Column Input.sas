/* raw data .txt file
5 obs, 4 vars

1234567890123456789012 ---- Columns indicator, not part of the data
Tim  M14510/21/1978
Sara  13009/20/1964
Mike M18011/23/1965
LauraF13011/06/1980
Sean M16704/07/2000

*/

data sdata_column;
   infile "/folders/myfolders/DATA_column.txt";
   input 	name	$  1-5
			Gender  $   6
			Weight	   7-9
			DOB     $ 10 - 19;
run;

proc print data = sdata_column;
run;

