# Table Code Breakdown

## Example for Observations Table:

```q
/ Observations Table
observations:([]
	obs_id:`long$();
	session_id:`int$();
	BJD_TDB:`float$();
	orbital_phase:`float$();
	flux:`float$();
	Uncertainty:`float$();
	model:`float$()
	airmass:`float$();
	amplitude`float$();
	offset:`float$();
	);
```

## Breakdown:

```q 
observations:([]...
```
* Assign name “observations” to table
### Table creation syntax:  
```q
tableName:([]
      col1:type$();  
      col2:type$(); ... 
      )
```

```q
 ([]...)
```
  * Creates an unkeyed table (plain table with no primary key)  
  * If you put column names inside the brackets like
    ```q
    (\[col1\] col2:...)
    ```
  * those become the keyed columns (like a primary key)

### Assignment   
  ```q
  tableName:(...)
  ```
    * `:` Assigns the result of the expression on the right to the variable on the left

### Type Casting 
```q
`type$()
```
  * Backtick symbol `\`int`  denoting the type name  
  * Cast operator `$` Casts what’s on its right to the type on the left  
  * List operator `()` → Empty list  
  * E.g.  Empty list of type integer
    ```q
    `int$()
    ```
  * See [Column Data Types doc](./Column_Data_Types.md) for more info on uses of data types

## Relational Structure (Foreign Keys):
  * The 5 tables form a hierarchy:

Sources (source\_id)  
→ obs\_sessions )sessions\_id, → source\_id)  
	→ observations (obs\_id, → session\_id)  
		→ centroids (--\> obs\_id) \[one-to-one\]  
		→ psf\_params (--\> obs\_id) \[one-to-one\]

* Sources \- master catalogue of stars/planets being studied  
* Obs\_sessions \- single night/session observing one source, with one telescope, one filter, one observer  
* Observations \- individual time-stamped flux measurements within a session (bulk data table)  
* Centroids \- pixel position of the star’s centre for each observation (one-to-one with observation)  
* Psf\_params \- Point Spread Function shape for each observation (one-to-one with observations)

- The foreign key relationships (source\_id, session\_id, obs\_id) are not enforced by the schema syntax. They’re relational by convention, enforced in queries.

* The Separator ;  
  * Separates expressions within a list of function bodies. Here it separates each column definition. The last column in each table has no trailing semicolon \- this is required q syntax (a trailing ; would add a null entry)

**Summary of what this code does**

* Creates 5 empty, typed, in-memory tables that form a relational schema for storing astronomical photometry data \- star/planet observations from telescopes, including the raw flux measurements, the sky positions, the pixel centroids, and the PSF shape parameters. The type choices (long for high-volume IDs, symbol for categorical text, float for scientific, date for temporal queries) are all deliberate performance decisions idiomatic to q/KDB+

