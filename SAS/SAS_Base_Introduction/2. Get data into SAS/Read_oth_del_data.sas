/*Import/Read data with delimiters other than Blanks or Commas*/

data sdata_colon;
   infile "/folders/myfolders/other_del_data.txt" dlm = ':'; 
/* 
   infile "/folders/myfolders/other_del_data.txt" delimiter = ':';  
   infile "/folders/myfolders/other_del_data.txt" dsd dlm = ':';    
*/
   input name $ Gender $ Age Weight;
run;
proc print data=sdata_colon;
run;

