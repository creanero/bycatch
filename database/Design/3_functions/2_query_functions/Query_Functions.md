**Query Functions**

**Query by time range:**  
*getObsByTimeRange:{\[start\_bjd;end\_bjd\]*  
	*obs:select from observations where BJD\_TDB within (start\_bjd;end\_bjd);*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};*

**Query by flux range:**  
*getObsByFlux:{\[min\_flux;max\_flux\]*  
	*obs:select from observations where flux within (min\_flux;max\_flux);*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};*

**Query by amplitude range:**  
*getObsByAmplitude:{\[min\_amp;max\_amp\]*  
	*obs:select from observations where amplitude within (min\_amp;max\_amp);*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};*

**Query by orbital phase range:**  
*getObsByPhase:{\[min\_phase;max\_phase\]*  
	*obs:select from observations where orbital\_phase within (min\_phase;max\_phase);*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};*

**Query by position (circular region on detector):**  
*getObsByPosition:{\[center\_x;center\_y;radius\]*  
	*cent:update dist:sqrt\[((x\_centroid-center\_x)xexp 2\) \+ (y\_centroid-center\_y)xexp 2\] from centroids;*  
	*cent\_filtered:select from cent where dist\<=radius;*  
	*observations lj \`obs\_id xkey cent\_filtered lj \`obs\_id xkey psf\_params*  
	*};*

* update add or modifies columns  
  * Adds new column dist (distance) computed from existing centroid coordinates  
  * Global centroids table is not modified \- update returns a new table stored in cent  
* xexp \- power operator:  
  * Raises left value to power of right, e.g. x xexp 2 is x2   
* Sqrt \= square root  
* All operations vectorised \- apply to every row of centroids table simultaneously, producing vector of distances  
* dest\<=radius \- filters only rows where distance is within search radius  
  * Join order: observations is left table and cent\_filtered is joined in rather than pre-filtered obs

**Get all data for a session:**  
*getSessionData:{\[sess\_id\]*  
	*obs:select from observations where session\_id=sess\_id;*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};*

* Uses \= instead of within \- exact match on foreign key session\_id rather than a range query

**Find transit events (flux drops below threshold):**  
*getTransitEvents:{\[start\_bjd;end\_bjd;max\_flux\]*  
	*obs:select from observations where BJD\_TDB within (start\_bjd;end\_bjd), flux\<max\_flux;*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};* 

* 2 conditions join by , (AND). Find observations in time range and below flux threshold

**Find high amplitude events in time range:**  
*getHighAmpEvents:{\[start\_bjd;end\_bjd;min\_amp\]*  
	*obs:select from observations where BJD within (start\_bjd;end\_bjd), amplitude\>min\_amp;*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};*

**Find high airmass observations (lower quality data):**  
*getHighAirmass:{\[min\_airmass\]*  
	*obs:select from observations where airmass\>min\_airmass;*  
	*obs lj \`obs\_id xkey centroids lj \`obs\_id xkey psf\_params*  
	*};*

* Simple filter: \> as threshold