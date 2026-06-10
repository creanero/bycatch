**Utility Functions**

**Display Statistics:**

*showStats:{\[\]*  
	*\-1 "=================================";*  
	*\-1 "Database Statistics:";*  
	*\-1 "  Sources: ",string count sources;*  
	*\-1 "  Sessions: ",string count obs\_sessions;*  
	*\-1 "  Observations: ",string count observations;*  
	*\-1 "  Centroids: ",string count centroids;*  
	*\-1 "  PSF Params: ",string count psf\_params;*  
	*\-1 "=================================";*  
	*};*

* *\-1* is the Print Operator  
  * Writes to standard output with a newline  
  * Equivalent to print() in Python  
* *showStats:{\[\]*  
  * No arguments, function takes nothing  
  * Reads global tables and prints  
* *\-1 “ Sources: “,string count sources;*  
  * Count sources → counts rows in sources table, returns integer  
* String → converts q value to string representation. E.g. string 42 gives “42”  
  * Can’t concatenate string and integer directly  
* “,” (comma) → join operator, concatenates 2 lists/strings together so “ Sources: “,”42” gives “Sources: 42”

**Save database to disc:**

*saveDB:{\[db\_path\]*  
	*(hsym \`$db\_path) set (sources;obs\_sessions;observations;centroids;psf\_params);*  
	*\-1 “Database saved to: “,db\_path;*  
	*};*

* *saveDB:{\[db\_path\]*  
  * Function taking 1 argument → string like “C:/kdb/db/test\_db” representing where to save the file  
* *(hsym \`$db\_path) set (sources;obs\_sessions;observations;centroids;psf\_params);*  
  * \\$db\_path\` → Casts string to symbol  
  * Hsym → Converts symbol into a file handle  
  * Set → Writes to disc  
* *(sources;obs\_sessions;observations;centroids;psf\_params)*  
  * *Plain list of tables, saving as list allows laodDB to retrieve them by position*  
    * *Index 0 for sources, index 1 for obs\_sessions, etc*

* *\-1 “Database saved to: “, db\_path*  
  * Prints confirmation to stout. , joins label string to path string

**Load database from disc:**

*loadDB:{\[db\_path\]*  
	*tables\_data:get hsym \`$db\_path;*  
	*\`sources set tables\_data 0;*  
	*\`obs\_sessions set tables\_data 1;*  
	*\`observations set tables\_data 2;*  
	*\`centroids set tables\_data 3;*  
	*\`psf\_params set tables\_data 4;*  
	*\-1 “Database loaded from: “,db\_path*  
	*};*

* *loadDB:{\[db\_path\]*  
  * 1 argument (db\_path) same string path passed to saveDB  
    * Must point to base filename  
* *tables\_data:get hsym \`db\_path;*  
  * \`$db\_path → String to symbol  
  * Hsym → converts symbol to file handle  
  * get → load operator → reads binary file from disc  
    * Stored in local variable tables\_data  
* *\`sources set tables\_data 0;* etc etc  
  * Symbol left of set assigns value to global variable with name rather than local one  
* *set*   
  * Symbol on left so assigns to variable rather than writing to disc  
* *tables\_data 0* etc etc  
  * Positional indexing  
    * Tables\_data is list (sources; etc;etc) as saved → index 0 \= sources, index 1 \= obs\_sessions etc etc.  
    * Order must match with saveDB  
* *\-1 “Database loaded from: “, db\_path*  
  * Prints confirmation to stout. , joins label string to path string

