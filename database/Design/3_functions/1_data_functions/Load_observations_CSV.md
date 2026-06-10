Data Functions: Load observations from CSV

1. *loadObservations:{\[csv\_path;session\_id\]*  
* Csv\_path → string filepath to CSV file on disk  
* Session\_id → Foreign key linking all these observations to their parent session (returned by a prior addSession call)  
2. *t:(“FFFFFFFFFFFFF”;enlist csv) 0: hsym \`$csv\_path;*  
* csv\_path → \`$csv\_path → hsym  
* \`$csv\_path  
  * \`$ casts string to symbol  
  * Csv\_path comes as string like “/data/obs.csv” and converts to \`/data/obs.csv  
* hsym \`$csv\_path  
  * Hsym converts symbol to file handle  
  * Results in file handle like :/data/obs.csv that can be opened and read  
* 0: file read operator → reads file handle on its right using format specs on its left  
* “FFFFFFFFFFFFF” → 13 F characters, 1 per column, putting every column as a float  
* Enlist csv → enlist wraps in list which signals to 0: that 1st row is header and has column names

3. *n:count t;*  
* Count returns number of rows in table t  
* Stored as n for reuse  
* Used to generate IDs and replicate the session\_id

4. *obs\_ids:NEXT\_OBS\_ID \+ til n;*  
   *NEXT\_OBS\_ID+:n;*  
* Til n → til generates a list of integers from 0 up to n (not incl. n)  
* NEXT\_OBS\_ID \+ til n  
  * Adding scalar to list is vectorised → applies to every element  
  * If NEXT\_OBS\_ID is 1 and n is 5, gives 1 2 3 4 5  
  * These become obs\_id values for all new rows  
* NEXT\_OBS\_ID+:n  
  * Increments counter by n → jumping forward by total number of rows just inserted so next load starts from correct ID

5. *\`observations upsert flip \`obs\_id\`session\_id\`BJD\_TDB\`orbital\_phase\`flux\`uncertainty\`model\`airmass\`amplitude\`offset\!(*   
* Backtick column names with no spaces → symbol list \- column names for new table being constructed  
* Keys\!values → creates dictionary mapping keys to values  
  * Symbol list of column names mapped to list of data vectors below it  
  * Results in dictionary of column → vector pair  
* Flip → transposes dictionary of column vectors into proper table  
  * E.g. flip \`col1\`col2\!(vector1;vector2) → table  
6. *n\#session\_id*  
* \# \= take operator → n\#session\_id repeats scalar session\_id value n times to create vector length n.   
  * Every observation row needs the same session\_id but column must be vector matching the length of all other columns  
7. *t\`BJD\_TDB*  
* Indexing table t with column name returns column as a vector.  
  * t\`BJD\_TDB gives all BJD\_TDB float values from the CSV as a list  
8. t\`$”Orbital Phase”  
* Added $”...” to account for space between words and avoid issues later on  
  * Converts string to symbol so q treats whole title as column name  
9. *\`centroids upsert flip \`obs\_id\`x\_centroid\`y\_centroid\!(*

*Obs\_ids;*  
*t\`x\_centroid;*  
*t\`y\_centroid*

* Same pattern  
* Same obs\_id vector (linking back to observations as foreign key) and the relevant columns pulled from t  
* Same for psf\_params as well  
10. *\`BJD\_TDB xasc \`observations;*

*\`obs\_id xasc \`centroids;*  
*\`obs\_id xasc \`psf\_params;*

* Xasc sorts table ascending by named column  
* Left argument is sort column (as symbol), right is table (as symbol, so global is modified in place)  
* Observations → sorted by \`BJD\_TDB \- chronological time order, natural order fro time-series analysis  
* Centroids & psf\_params sorted by obs\_id \- keeping aligned with observations  
11. *n* → Return count  
*  Returns the count of rows inserted \- confirm how many observations were loaded

