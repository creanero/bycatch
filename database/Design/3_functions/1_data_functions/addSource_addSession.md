**Data Functions: addSource & addSession**

**Function Syntax**

*functionName:{\[arg1;arg2;arg3\]*  
	*/ body*  
	*}* 

* { and } → curly brackets define lambda (function) in q  
* \[arg1;arg2;arg3\] → argument list, inside square brackets, separated by ; (same separator rule as columns in tables)  
* The function body is everything between the brackets and the closing }  
* The last expression in the body is always the return value \- q has no explicit return keyword

**Add a new source / Add new observation session**

Add a new source:

1. *addSource:{\[name;ra\_val;dec\_val;note\]*  
* Name \- star/planet name (will become a symbol)  
* Ra\_val \- Right ascension float value  
* Dec\_val \- Declination float value  
* Note \- notes symbol  
2. *source\_id:NEXT\_SOURCE\_ID;*  
* Captures current value of global counter into local variable source\_id  
* This is ID that will be assigned to new row  
* Saved first before counter is incremented, so ID is frozen at this moment  
3. *NEXT\_SOURCE\_ID+:1;*  
* Compound assignment operator   
* \+: means  increment in place  
* Equivalent to NEXT\_SOURCE\_ID: NEXT\_SOURCE\_ID \+ 1  
* Counter is bumped so next call to addSource gets next ID  
4. *\`sources upsert (source\_id;name;ra\_val;dec\_val;note);*  
* \`sources → backtick symbol referencing table by name  
* Using backtick (symbol reference) rather than bare name sources tells q to modify the global table in place, rather than a local copy  
* Upsert → q’s combined insert/update operator  
* For unkeyed table, like sources, always appends a new row  
* On a keyed table it would update if the key exists, insert if not  
5. *(sources\_id;name;ra\_val;dec\_val;note)*  
* Simple list of values, one per column, in same order columns were defined in the table schema. q matches them positionally.  
6. *Sources\_id*  
* The return value → last expression with no ; terminator  
* Function returns newly assigned sources\_id so the caller knows what ID was created  
* Important for then calling addSession, which needs to pass in a source\_id as a foreign key

**Add a new observation session:**

1. *addSession:{\[src\_id;tel;filt;obs;cond;sess\_date;note\]*  
* Src\_id → foreign key to sources.sources\_id (returned by prior addSource call)  
* tel → telescope name symbol  
* Filt → filter symbol (e.g. \`R, \`V)  
* Obs → observer name symbol  
* Cond → conditions symbol (e.g. \`good)  
* Sess\_date → date value  
* Note → notes symbol

* 2 functions establish intended call order that mirrors table hierarchy  
* *sid: addSource\[\`WASP-39b; 330.67; \-3.44; \`transit-target\]; / returns e.g. 1 eid: addSession\[sid; \`Kepler; \`R; \`Smith; \`good; 2024.06.01; \`clear\]; / uses sid=1*   
* The returned ID from addSource feeds directly into addSession as src\_id, maintaining the foreign key relationship manually \- q doesn’t enforce this automatically