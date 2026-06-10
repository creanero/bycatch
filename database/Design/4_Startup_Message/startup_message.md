**Start-up Message**

*\-1 "\\n========================================";*  
*\-1 "Astronomical Bycatch Database Loaded";*  
*\-1 "========================================";*  
*\-1 "\\nQuick Start (HAT-P-32b example):";*  
*\-1 "  1\. source\_id:addSource\[\`$\\”HAT-P-32b\\”; 0n; 0n; \`$\\"Hot Jupiter\\"\]";*  
*\-1 "  2\. session\_id:addSession\[source\_id; \`Kepler; \`R; \`Pipeline; \`good; 2018.01.01; \`\\"\\"\]";*  
*\-1 "  3\. n:loadObservations\[\\"HAT-P-32b.csv\\"; session\_id\]";*  
*\-1 "  4\. showStats\[\]";*  
*\-1 "\\nQuery Functions:";*  
*\-1 "  getObsByTimeRange\[start\_bjd; end\_bjd\]     / Query by BJD\_TDB";*  
*\-1 "  getObsByFlux\[min\_flux; max\_flux\]          / Query by normalized flux";*  
*\-1 "  getObsByAmplitude\[min\_amp; max\_amp\]       / Query by amplitude";*  
*\-1 "  getObsByPhase\[min\_phase; max\_phase\]       / Query by orbital phase";*  
*\-1 "  getObsByPosition\[x; y; radius\]            / Query by detector position";*  
*\-1 "  getSessionData\[session\_id\]                / Get full session data";*  
*\-1 "  getTransitEvents\[start; end; max\_flux\]    / Find transit dips";*  
*\-1 "  getHighAirmass\[min\_airmass\]               / Find high airmass obs";*  
*\-1 "\\nDatabase Management:";*  
*\-1 "  saveDB\[\\"path/to/save\\"\]    / Save to disc";*  
*\-1 "  loadDB\[\\"path/to/load\\"\]    / Load from disc";*  
*\-1 "========================================\\n";*

* Executes immediately as soon as script is loaded  
* \\n is newline escape sequence  
  * Inserts blank line before \=== separator  
* String literal to show user example code to type into database  
  * \` (backtick) inside string \-\> literal characters being printed  
  * Showing user what to type  
* \\” → Escaped quotes. E.g \`\\”Hot Jupiter\\”  
  * Because whole thing wrapped in double quotes, any double quote that needs to appear inside string must be escaped with backlash  
  * So \`\\”Hot Jupiter\\” prints as \`”Hot Jupiter” \- showing user correct syntax for multi-word symbol  
* 0n (Null float)  
  * Valid placeholder for values if not known/being used yet

* Quick Start example deliberately demonstrates the correct call order and shows how the return values chain together:  
  * addSource → returns source)id → passed into addSession  
  * addSession → returns session\_id → passed into loadObservations

***Example file path for save and load DB:***

* saveDB\[“C:/kdb/db/test2\_db”\]  
* loadDB\[“C:/kdb/db/test2\_db”\]