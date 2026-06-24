//// Astronomical Bycatch Database ////
/ Ensure the current directory contains this file before loading
/ Load this file with: \l targetStarOnly.q


/ Section to create tables of for input data with specified columns and data types
/// Create Tables
/ Stores info about astronomical sources (stars)
// Sources Table
sources:([]				/ Create unkeyed table for sources data - column names inside the brackets give keyed columns
	source_id:`int$();		/ Unique integer ID for each star/planet system (Use star designation as ID) (int allows ~2 billion sources)
	source_name:`symbol$();		/ Star/planet name as symbol (e.g. HAT
	ra:`float$();			/ Right ascension in degrees (float for decimal precision in sky coordinates)
	dec:`float$();			/ Declination in degrees (float for decimal precision in sky coordinates)
	notes:`symbol$()		/ Additional info as symbol (efficient storage for repeated text values) [No ";" after last column]
	);				/ ")" Closes table definition and ";" terminates the whole statement, seperating it from any following script

/ Stores metadata about each observing session
/ One session = one period of telescope observation
/ Links to sources via source_id
	/ One source can have many sessions
// Observation Sessions Table
obs_sessions:([]			/ Create unkeyed table for observation sessions data - column names inside the brackets give keyed columns
	session_id:`int$();		/ Unique integer ID for each observing session (int sufficient for typical session counts)
	source_id:`int$();		/ Foreign key to sources table (int matches source_id type for joins)
	telescope:`symbol$();		/ Telescope name as symbol (efficient for repeated values like Kepler, Hubble)
	filter:`symbol$();		/ Optical filter used as symbol (e.g. R, V, I - limited distinct values)
	observer:`symbol$();		/ Observer name/ID as symbol (efficient for repeated names)
	conditions:`symbol$();		/ Observing conditions as symbol (e.g. good, poor - categorical data)
	session_date:`date$();		/ Calendar date of observation (date type for temporal queries and sorting)
	notes:`symbol$()		/ Additional info as symbol (efficient storage for repeated text values)
	);				/ ")" Closes table definition and ";" terminates the whole statement, seperating it from any following script

/ Stores the photometric measurements
/ Links to obs_sessions via session_id
	/ One session has many observations
// Observations Table
observations:([]			/ Create unkeyed table for observations data - column names inside the brackets give keyed columns
	obs_id:`long$();		/ Unique observations ID as long (64-bit allows billions of measurements)
	session_id:`int$();		/ Foreign key to obs_sessions table (int matches session_id type for joins)
	BJD_TDB:`float$();		/ Barycentric Julian Date in TDB (Barycentric Dynamical Time) timescale (float for fractional phase calculations)
	orbital_phase:`float$();	/ Planet's position in orbit 0-1 (float for fractional phase calculations)
	flux:`float$();			/ Normalized stellar flux (float for precise photometric ratios)
	uncertainty:`float$();		/ Photometric uncertainty/error (float for statistical analysis and error propagation)
	model:`float$();		/ Theoretical model flux value (float to match observed flux precision)
	airmass:`float$();		/ Atmospheric path length multiplier (float for extinxtion corrections)
	amplitude:`float$();		/ Stellar brightness amplitude (float for photometric measurements)
	offset:`float$()		/ Calibration offset applied to data (float for additive corrections)
	);				/ ")" Closes table definition and ";" terminates the whole statement, seperating it from any following script

/ Stores star's poistion on detector
/ Links 1:1 with observations via obs_id
/ Used for quality control and aperture photometry
// Centroids Table
centroids:([]				/ Create unkeyed table for centroids data - column names inside the brackets give keyed columns
	obs_id:`long$();		/ Foreign key to observations table (long matches obs_id for one-to-one linking)
	x_centroid:`float$();		/ Horizontal pixel position of star centre (float for sub-pixel centroiding precision)
	y_centroid:`float$()		/ Vertical pixel position of star centre (float for sub-pixel centroiding precision)
	);				/ ")" Closes table definition and ";" terminates the whole statement, seperating it from any following script

/ Stores Point Spread Function data
/ Links 1:1 with observations via obs_id
/ Describes how blurred/spread out teh star image appears
// PSF Parameters Table
psf_params:([]				/ Create unkeyed table for PSF parameters data - column names inside the brackets give keyed columns
	obs_id:`long$();		/ Foreign key to observations table (long matches obs_id for one-to-one linking)
	sigma_x:`float$();		/ Gaussian PSF width in x-direction in pixels (float for sub-pixel fitting precision)
	sigma_y:`float$();		/ Gaussian PSF width in y-direction in pixels (float for sub-pixel fitting precision)
	rotation:`float$()		/ PSF rotation angle in radians (float for precise angular measurements)
	);				/ ")" Closes table definition and ";" terminates the whole statement, seperating it from any following script


/// ID Counter Initialisation
/ Global counters for auto-generating unique IDs
/ These increment each time a new record is created
// ID Counters
/ *** Note: NEXT_SOURCE_ID not used since source_id is symbol type (uses star name directly) ***
NEXT_SOURCE_ID:1;	/ Auto-increment counter for source_id - starts at 1, increments with each new source inserted
NEXT_SESSION_ID:1;	/ Auto-increment counter for obs_session_id - starts at 1, increments with each new session inserted
NEXT_OBS_ID:1;		/ Auto-increment counter for session_id - starts at 1, increments with each new session inserted


/// Data Functions

// Add a new source
addSource:{[name;ra_val;dec_val;note]				/ Function taking 4 args
	source_id:NEXT_SOURCE_ID;				/ Capture current counter value as this row's ID before incrementing
	NEXT_SOURCE_ID+:1;					/ Increment global counter in place (+: is compound add-assign) so next call gets next ID
	`sources upsert (source_id;name;ra_val;dec_val;note);	/ Upsert positional list into global sources table (backtick = modify global in place)
	source_id						/ Return new ID to caller (no semicolon = implicit return of last expression)
	};

// Add a new observation session
addSession:{[src_id;tel;filt;obs;cond;sess_date;note]					/ Function taking 7 args 
	session_id:NEXT_SESSION_ID;							/ Capture current counter value as this row's ID before incrementing 
	NEXT_SESSION_ID+:1;								/ Increment global counter in place so next call gets ID
	`obs_sessions upsert (session_id;src_id;tel;filt;obs;cond;sess_date;note);	/ Upsert positional list into global obs_sessions table
	session_id									/ Return the new session ID to caller for use as foreign key in addObservation calls
	};

// Load observations from CSV
/ CSV must have columns: BJD_TDB, Orbital Phase, Flux, Uncertainity, Model, Airmass, x_centroid, y_centroid, amplitude, sigma_x, sigma_y, rotation, offset
/ First line should be the header
loadObservations:{[csv_path;session_id]
	/ Load CSV file with header
	/ 13 float columns
	/ enlist csv = use comma delimiter & treat 1st row as column headers
	/ 0: = file read operator
	/ hsym `$csv_path = cast string path to symbol then to file handle to open
	t:("FFFFFFFFFFFFF";enlist csv) 0: hsym `$csv_path;
    
	/ Number of observations  
	n:count t;				/ Count rows in loaded table -> stored for reuse in ID generation and session_id replication
    
	/ Generate unique obs_ids
	/ Generat n sequential unique obs_ids starting from current counter value
	/ til n -> produces 0,1,2...n-1 and adding NEXT_OBS_ID shifts range to NEXT_OBS_ID, NEXT_OBS_ID+1...
	obs_ids:NEXT_OBS_ID + til n;
	NEXT_OBS_ID+:n;				/ Advance global counter by n - to account for all rows inserted
    
	/ Build & Insert observations table
	/ Column names match CSV header exactly
	/ flip converts dictionary (symbol_list!vector_list) into proper table
	/ n#session_id replicates scalar session_id into a vector of length n (one value per row)
	/ `t colname extracts named column from loaded CSV table as a vector
	`observations upsert flip `obs_id`session_id`BJD_TDB`orbital_phase`flux`uncertainty`model`airmass`amplitude`offset!(
		obs_ids;			/ Unique obs_ids generated to fill column
		n#session_id;			/ Session foreign key replicated n times to fill column
		t`BJD_TDB;			/ Barycentric julian date vetcor from CSV
		t`orbital_phase;		/ Orbital phase vector from CSV
		t`flux;				/ Normalised flux vector from CSV
		t`uncertainty;			/ Photometric uncertainty
		t`model;			/ Theoretical model flux vector from CSV
		t`airmass;			/ Airmass vector from CSV
		t`amplitude;			/ Brightness amplitude vector from CSV
		t`offset			/ Calibration offset vector from CSV
	);

	/ Build & Insert into centroids table - using same obs_ids to maintain foreign key link to observations
 	`centroids upsert flip `obs_id`x_centroid`y_centroid!(
 		obs_ids;			/ same obs_ids -> one-to-one foreign key link to observations
 		t`x_centroid;			/ Horizontal pixel centroid vector from CSV
		t`y_centroid			/ Vertical pixel centroid vector from CSV
 	);
    
	/ Build & Insert into psf_params table using same obs_ids to maintain foreign key link to observations
	`psf_params upsert flip `obs_id`sigma_x`sigma_y`rotation!(
		obs_ids;			/ same obs_ids -> one-to-one foreign key link to observations
		t`sigma_x;			/ PSF x-width vector from CSV
		t`sigma_y;			/ PSF y-width vector from CSV
		t`rotation			/ PSF rotation angle vector from CSV
	);
    
	/ Sort tables
	`BJD_TDB xasc `observations;		/ Sort observations chronologically by barycentric julian date
	`obs_id xasc `centroids;		/ Sort centroids by obs_id to align with observations
	`obs_id xasc `psf_params;		/ Sort psf_params by obs_id to align with observations
    
	/ Return count
	n					/ Return number of rows inserted as confirmation
	};


/// Query Functions

// Query by time range
/ Return all observations (with centroids and PSF data) within a BJD time range
getObsByTimeRange:{[start_bjd;end_bjd]
	obs:select from observations where BJD_TDB within (start_bjd;end_bjd);	/ within = inclusive range filter, returns all 
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params		/ Chain lj left joins - xkey creates temporary keyed table for join, evaluated right to left
	};

// Query by flux range
/ Return all observations within a normalised flux range
getObsByFlux:{[min_flux;max_flux]
	obs:select from observations where flux within (min_flux;max_flux);
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Query by amplitude range
/ Return all observations within a brightness amplitude range
getObsByAmplitude:{[min_amp;max_amp]
	obs:select from observations where amplitude within (min_amp;max_amp);
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Query by orbital phase range
/ Return all observations within an orbital phase range (0-1)
getObsByPhase:{[min_phase;max_phase]
	obs:select from observations where orbital_phase within (min_phase;max_phase);
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Query by position (circular region on detector)
/ Return all observations where the star centroid falls within a cicular region on the detector
getObsByPosition:{[center_x;center_y;radius]
	/ update adds computed column dist to a new table - doesn't modify global centroids
	/ xexp = power operator, sqrt = square root - implements 2D Pythagorean distance formula across all rows vectorised
	cent:update dist:sqrt[((x_centroid-center_x)xexp 2) + (y_centroid-center_y)xexp 2] from centroids;
	cent_filtered:select from cent where dist<=radius;				/ Filter to rows within circular radius
	observations lj `obs_id xkey cent_filtered lj `obs_id xkey psf_params		/ Join filtered centorids and psf onto full observations table
	};

// Get all data for a session
/ Return all observations belonging to s specific session (exact foreign key match)
getSessionData:{[sess_id]
	obs:select from observations where session_id=sess_id;		/ = exact match rather than within range
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Find transit events (flux drops below threshold)
/ Return observations in time range where flux drops below threshol - identifies planetary transit events
getTransitEvents:{[start_bjd;end_bjd;max_flux]
	obs:select from observations where BJD_TDB within (start_bjd;end_bjd), flux<max_flux;	/ comma = AND -> both conditions must be true
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Find high amplitude events in time range
/ Return high brightness amplitude events in a time range
getHighAmpEvents:{[start_bjd;end_bjd;min_amp]
	obs:select from observations where BJD_TDB within (start_bjd;end_bjd), amplitude>min_amp;
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Find high airmass observations (lower quality data)
/ Return observations with high airmass indicating lower quality data (longer atmospheric path)
getHighAirmass:{[min_airmass]
	obs:select from observations where airmass>min_airmass;		/ Single threshold filter, no range needed
	obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};


/// Utility Functions

// Display statistics
/ Display row counts for all 5 tables as a formatted summary
showStats:{[]							/ No arguments - read global tables directly
	-1 "=================================";			/ -1 = Write to stdout with newline
	-1 "Database Statistics:";
	-1 "  Sources: ",string count sources;			/ Count rows -> convert int to string -> join to label + print
	-1 "  Sessions: ",string count obs_sessions;
	-1 "  Observations: ",string count observations;
	-1 "  Centroids: ",string count centroids;
	-1 "  PSF Params: ",string count psf_params;
	-1 "=================================";
	};

// Save database to disc
/ Serialise all 5 tables to single q binary file at db_path
saveDB:{[db_path]
	/ hsym `$db_path: cast string path to symbol (`$) then to file handle (hsym) so set knows to write to disc
	/ set: writes data to disc as q native binary fomrat(splits into multiple files automatically for large data)
	/ (sources;obs_sessions;...): saves tables as a positional list so loadDB can retrieve by index (0=sources, 1_obs_sessions etc.)
	(hsym `$db_path) set (sources;obs_sessions;observations;centroids;psf_params);	
	-1 "Database saved to: ",db_path;	/confirm save to stdout - , joins label string to path string
	};

// Load database from disc
/ Load previously saved q binary file back into meory, restoring all tables as globals
loadDB:{[db_path]
	/ get: reads q binary from disc - mirro of set used in saveDb
	/ hsym `$db_path: convert string path to file handle so get knows to read from disc
	/ Result is positional list as saved
	tables_data:get hsym `$db_path;
	`sources set tables_data 0;		/ set with symbol on left assigns to global variabel; index 0 = sources table
	`obs_sessions set tables_data 1;	/ same for rest --> Order must match saveDB exactly
	`observations set tables_data 2;
	`centroids set tables_data 3;
	`psf_params set tables_data 4;
	-1 "Database loaded from: ",db_path;	/ Confirm load to stout by joining path string to label
	};


/// Start-up Message
/ Top-level code -> executes immediately on script load

-1 "\n========================================";	/ \n = newlinw escape - blank line before banner
-1 "Astronomical Bycatch Database Loaded";		/ Script identity banner
-1 "========================================";

-1 "\nQuick Start (HAT-P-32b example):";							/ \n adds blank line before each section heading
-1 "  1. source_id:addSource[`$\"HAT-P-32b\"; 0n; 0n; `$\"Hot Jupiter\"]";			/ 0n = null float placeholder for unknown coords
-1 "  2. session_id:addSession[source_id; `Kepler; `R; `Pipeline; `good; 2018.01.01; `none]";	/ Shows how returned source_id feeds into addSession as foreign key
-1 "  3. n:loadObservations[\"HAT-P-32b.csv\"; session_id]";					/ Shows how returned session_id feeds into observations as foreign key
-1 "  4. showStats[]";										/ Shows utility function call with no arguments

-1 "\nQuery Functions:";
/ Documentation strings - / insidde is printed literally, not a comment 
-1 "  getObsByTimeRange[start_bjd; end_bjd]     / Query by BJD_TDB";
-1 "  getObsByFlux[min_flux; max_flux]          / Query by normalized flux";
-1 "  getObsByAmplitude[min_amp; max_amp]       / Query by amplitude";
-1 "  getObsByPhase[min_phase; max_phase]       / Query by orbital phase";
-1 "  getObsByPosition[x; y; radius]            / Query by detector position";
-1 "  getSessionData[session_id]                / Get full session data";
-1 "  getTransitEvents[start; end; max_flux]    / Find transit dips";
-1 "  getHighAirmass[min_airmass]               / Find high airmass obs";

-1 "\nDatabase Management:";
-1 "  saveDB[\"path/to/save\"]    / Save to disc";	/ \" escapes double quote inside outer string delimiters
-1 "  loadDB[\"path/to/load\"]    / Load from disc";
-1 "========================================\n";	/ Trailing \n adds blank line after banner for clean terminal prompt
