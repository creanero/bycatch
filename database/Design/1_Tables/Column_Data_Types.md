**Column Data Types**

* int → Integer  
  * Used for IDs and Foreign Keys  
  * \~ \+/- 2.1 billion  
* long   
  * Used for ‘obs\_id’ where billions of rows are expected  
  * \~ \+/- 9.2 quintillion  
* float   
  * Used for coordinates, measurements, scientific values  
  * Decimal precision  
* symbol   
  * Interned string  
  * Stored as a pointer to a single copy \- efficient for repeated text like “good”, “Kepler”  
* date   
  * Days since 2000.01.01  
  * Efficient for temporal queries 

* Why symbol for text?  
  * Symbols are interned  
    * If 1,000 rows all have “good” as their conditions, they all point to the same memory location  
    * This is much more efficient than storing the string “good” 1,000 times  
    * Categorical/Repeated values like telescope names, observer names, filter types, and conditions all use the data type of symbol

* Why long for obs\_id?  
  * Observations can be in the billions for a telescope producing continuous photometric data  
    * Int tops out at \~2.1 billion  
    * Long gives \~9.2 quintillion  
    * Obs\_id foreign key in centroids and psf\_params also uses long to match types exactly → q is strict about type consistency in joins

* Why float for scientific measurements?  
  * BJD\_TDB, flux, RA, DEC, etc all need sub-decimal precision  
  * q’s float give \~15-16 significant figures → sufficient for astronomical coordinates and photometric rations