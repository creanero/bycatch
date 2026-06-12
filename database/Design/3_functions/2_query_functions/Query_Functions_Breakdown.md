**Query Functions Breakdown**

* Query functions filter and join data across observation related tables  
  * 2 main concepts are *select* and *lj* join

* select \- query syntax:  
  * Structure: select \<columns\> from \<table\> where \<conditions\>  
  * If you omit column list, it returns all columns (equivalent to SELECT \* in SQL  
  * The result is always a new table

* within \- range operator:  
  * Within tests whether a value falls inclusively between two bounds  
  * (start\_bjd;end\_bjd) is a two-element list defining the range  
  * \[Equivalent in SQL is BETWEEN start AND end\]  
  * Works on any comparable type \- floats, dates, integers

* lj \- left join:  
  * Left join operator \- keeps all rows from left table and joins matching rows from the right table on a shared key column  
  * Non-matching rows get nulls for the right table’s columns \[same as SQL LEFT JOIN\]  
  * Left table (obs) already filtered, Right table must be keyed before joining

* xkey \- keying a table:  
  * Xkey sets the key column(s) of a table, returning a keyed copy  
  * Left argument is column name(s) as a symbol, right is the table  
  * lj requires right hand table to be keyed so it knows which column to join on  
  * E.g. *\`obs\_id xkey centroids* → Not modifying global centroids table 0 it’s creating a temporary keyed version just for this join operation

* Chained joins:

E.g. *obs lj \`obs\_id xkey centroids lj \`obs)id xkey psf\_params*

* Evaluated right to left  
  * 1\. \`obs\_id xkey psf\_parasm → keys psf\_params on obs\_id  
  * 2\. Centroids lj the result → joins psf\_parasm columns onto centroids  
  * 3\. \`obs\_id xkey that result → keys the combined centrodis+psf table  
  * 4\. obs lj that result → joins everything onto the filtered observations  
* Results in wide table combining all columns from observations, centroids, and psf\_parasm for the filtered rows