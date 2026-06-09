**ID Counters Code Breakdown**

**Full Code:**

*NEXT\_SOURCE\_ID:1;*	  
*NEXT\_SESSION\_ID:1;*	  
*NEXT\_OBS\_ID:1;*

**Example for NEXT\_SESSION\_ID:**

* NEXT\_SESSION\_ID (Variable name)  
  * All caps is naming convention to signal global constant/counter, distinct from table names and column names  
* “:” is the assignment operator in q  
* “1” an integer literal  
  * Counters start at 1; 1st record inserted into observations will get obs\_id of 1, the 2nd will get 2, etc  
* “;” is the statement terminator

**Note:**

* Counters are, currently, in-memory global variables. If q is restarted, they reset to 1\.   
* In production system, either persist to disk or derive next ID dynamically from the table itself (e.g. 1 \+ count sources)
