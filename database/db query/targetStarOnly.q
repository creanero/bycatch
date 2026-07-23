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

// Target Star Photometry Table 
targetStarOnly:([]
    id:`int$();
    time:`float$();
    T_flux:`float$();
    diff_T_flux:`float$();
    C1:`float$();
    C2:`float$();
    comp_star_check:`float$();
    reference_flux:`float$();
    t_error:`float$();
    diff_error:`float$();
    norm_t:`float$();
    oot_mask:`boolean$();
    norm_flux:`float$();
    norm_flux_error:`float$()
    );

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
/ Observations Table
observations:([]			/ Create unkeyed table for observations data - column names inside the brackets give keyed columns
	obs_id:`long$();		/ Unique observations ID as long (64-bit allows billions of measurements)
	session_id:`int$();		/ Foreign key to obs_sessions table (int matches session_id type for joins)
	BJD_TDB:`float$();		/ Barycentric Julian Date in TDB (Barycentric Dynamical Time) timescale (float for fractional phase calculations)
	flux:`float$();			/ Stellar flux measurement
	uncertainty:`float$();		/ Photometric uncertainty/error
	norm_t:`float$();		/ Normalized time coordinate used for phase folding
	oot_mask:`boolean$();		/ Out-of-transit mask flag
	norm_flux:`float$();		/ Normalized flux after baseline correction
	norm_flux_error:`float$();	/ Error on the normalized flux
	orbital_phase:`float$();	/ Orbital phase value when available
	model:`float$();		/ Model fit value when available
	airmass:`float$();		/ Airmass when available
	amplitude:`float$();	/ Amplitude parameter when available
	offset:`float$()		/ Offset parameter when available
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

// Records each CSV load event for traceability
target_load_log:([]
	target_name:`symbol$();
	csv_path:`symbol$();
	session_id:`int$();
	rows_loaded:`int$();
	loaded_at:`timestamp$()
	);


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

// Infer a target symbol from a CSV path (e.g. HATP-32-b.csv -> `HATP-32-B)
inferTargetFromCSV:{[csv_path]
	path_norm:ssr[csv_path;"\\";"/"];
	fname:last "/" vs path_norm;
	stem:first "." vs fname;
	norm:upper ssr[ssr[stem;"_";"-"];" ";""];
	`$norm
	};

// Reuse source rows when they already exist, otherwise create a new source.
getOrCreateSource:{[target_name]
	existing:select source_id from sources where source_name=target_name;
	$[0<count existing; first existing`source_id; addSource[target_name;0n;0n;`$"Auto-loaded target"]
	]
	};

// Create a standard default session for automated CSV ingestion.
createDefaultSession:{[source_id]
	addSession[source_id;`unknown;`unknown;`pipeline;`unknown;.z.D;`auto]
	};

// Ingest photometry or EXOTIC-style output rows into the database schema
/ Accept either a q table or a CSV file path.
/ This loader maps whichever columns are present to the shared schema.
loadObservations:{[input;session_id]
	/ If a string path is supplied, read the CSV into a q table.
	/ Otherwise assume the input is already a q table.
	t:$[10h=type input; 0!(("IFFFFFFBFF";enlist csv) 0: hsym `$input); 0!input];
	n:count t;

	/ Resolve the available column names from either the photometry pipeline or EXOTIC-style inputs.
	resolveCol:{[tbl;candidates]
		cols_t:cols tbl;
		check:{[cols_t;x] x in cols_t};
		present:candidates where check[cols_t] each candidates;
		$[0<count present; first present; `]
	};
	timeCol:resolveCol[t;(`time;`BJD_TDB)];
	fluxCol:resolveCol[t;(`T_flux;`Flux;`flux)];
	uncertCol:resolveCol[t;(`uncertainty;`Uncertainty)];
	normTCol:resolveCol[t;enlist `norm_t];
	ootMaskCol:resolveCol[t;enlist `oot_mask];
	normFluxCol:resolveCol[t;enlist `norm_flux];
	normFluxErrCol:resolveCol[t;enlist `norm_flux_error];
	phaseCol:resolveCol[t;(`orbital_phase;`$"Orbital Phase")];
	modelCol:resolveCol[t;(`model;`Model)];
	airmassCol:resolveCol[t;(`airmass;`Airmass)];
	amplitudeCol:resolveCol[t;enlist `amplitude];
	offsetCol:resolveCol[t;enlist `offset];
	xCol:resolveCol[t;enlist `x_centroid];
	yCol:resolveCol[t;enlist `y_centroid];
	psfXCol:resolveCol[t;enlist `sigma_x];
	psfYCol:resolveCol[t;enlist `sigma_y];
	rotCol:resolveCol[t;enlist `rotation];

	/ Explicitly cast the incoming values to the target schema types.
	obs_ids: NEXT_OBS_ID + til n;
	NEXT_OBS_ID+: n;
	obs_session_ids: n#session_id;
	bjd_vals: $[timeCol~`; n#0n; `float$ t timeCol];
	flux_vals: $[fluxCol~`; n#0n; `float$ t fluxCol];
	uncertainty_vals: $[uncertCol~`; n#0n; `float$ t uncertCol];
	norm_t_vals: $[normTCol~`; n#0n; `float$ t normTCol];
	oot_mask_vals: $[ootMaskCol~`; n#0n; `boolean$ t ootMaskCol];
	norm_flux_vals: $[normFluxCol~`; n#0n; `float$ t normFluxCol];
	norm_flux_error_vals: $[normFluxErrCol~`; n#0n; `float$ t normFluxErrCol];
	orbital_phase_vals: $[phaseCol~`; n#0n; `float$ t phaseCol];
	model_vals: $[modelCol~`; n#0n; `float$ t modelCol];
	airmass_vals: $[airmassCol~`; n#0n; `float$ t airmassCol];
	amplitude_vals: $[amplitudeCol~`; n#0n; `float$ t amplitudeCol];
	offset_vals: $[offsetCol~`; n#0n; `float$ t offsetCol];
	x_vals: $[xCol~`; n#0n; `float$ t xCol];
	y_vals: $[yCol~`; n#0n; `float$ t yCol];
	sigma_x_vals: $[psfXCol~`; n#0n; `float$ t psfXCol];
	sigma_y_vals: $[psfYCol~`; n#0n; `float$ t psfYCol];
	rotation_vals: $[rotCol~`; n#0n; `float$ t rotCol];

	/ Insert into observations table.
	`observations upsert flip `obs_id`session_id`BJD_TDB`flux`uncertainty`norm_t`oot_mask`norm_flux`norm_flux_error`orbital_phase`model`airmass`amplitude`offset!(
		obs_ids;
		obs_session_ids;
		bjd_vals;
		flux_vals;
		uncertainty_vals;
		norm_t_vals;
		oot_mask_vals;
		norm_flux_vals;
		norm_flux_error_vals;
		orbital_phase_vals;
		model_vals;
		airmass_vals;
		amplitude_vals;
		offset_vals
	);

	/ Insert into centroids table.
	`centroids upsert flip `obs_id`x_centroid`y_centroid!(
		obs_ids;
		x_vals;
		y_vals
	);

	/ Insert PSF parameters when available; otherwise use nulls.
	`psf_params upsert flip `obs_id`sigma_x`sigma_y`rotation!(
		obs_ids;
		sigma_x_vals;
		sigma_y_vals;
		rotation_vals
	);

	/ Sort tables for consistent retrieval order.
	`BJD_TDB xasc `observations;
	`obs_id xasc `centroids;
	`obs_id xasc `psf_params;

	/ Return the number of rows ingested.
	n
	};

// Load CSV with explicit target name (preferred when file names are generic)
loadTargetCSVAs:{[csv_path;target_name]
	target_sym:`$target_name;
	source_id:getOrCreateSource[target_sym];
	session_id:createDefaultSession[source_id];
	n:loadObservations[csv_path; session_id];
	`target_load_log upsert (target_sym;`$csv_path;session_id;`int$n;.z.P);
	show n;
	showStats[];
	n
	};

// One-line helper to infer target from file name, then load into DB
loadTargetCSV:{[csv_path]
	target_name:inferTargetFromCSV[csv_path];
	loadTargetCSVAs[csv_path;target_name]
	};

printFITSHint:{[]
	-1 "FITS input detected.";
	-1 "Run photometry pipeline first to generate CSV, then call sortTargetInput on that CSV.";
	-1 "Example: .\\.venv\\Scripts\\python db_pipeline.py target --csv my_output.csv";
	0N
	};

printUnsupportedInputHint:{[]
	-1 "Unsupported input type. Provide a CSV path or a FITS/FIT/FTS path.";
	0N
	};

// Generic entrypoint for file-based inputs.
// CSV is ingested directly; FITS is routed to the photometry pipeline first.
sortTargetInput:{[input_path;target_name]
	path_norm:ssr[input_path;"\\";"/"];
	parts:"." vs path_norm;
	ext:lower last parts;
	$[ext~"csv";
		loadTargetCSVAs[input_path;target_name];
		(ext~"fits") or (ext~"fit") or (ext~"fts");
		printFITSHint[];
		printUnsupportedInputHint[]
		]
	};


/// Query Functions

// Query by time range
/ Return all observations (with centroids and PSF data) within a BJD time range
getObsByTimeRange:{[start_bjd;end_bjd]
	obs:select from observations where BJD_TDB within (start_bjd;end_bjd);	/ within = inclusive range filter, returns all 
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params		/ Chain lj left joins - xkey creates temporary keyed table for join, evaluated right to left
	};

// Query all observations for a specific target (string or symbol)
getObsByTarget:{[target_name]
	target_sym:`$target_name;
	target_source_ids:exec source_id from sources where source_name=target_sym;
	target_session_ids:exec session_id from obs_sessions where source_id in target_source_ids;
	obs:select from observations where session_id in target_session_ids;
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Query by flux range
/ Return all observations within a normalised flux range
getObsByFlux:{[min_flux;max_flux]
	obs:select from observations where flux within (min_flux;max_flux);
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Query by amplitude range
/ Return all observations within a brightness amplitude range
getObsByAmplitude:{[min_amp;max_amp]
	obs:select from observations where amplitude within (min_amp;max_amp);
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Query by orbital phase range
/ Return all observations within an orbital phase range (0-1)
getObsByPhase:{[min_phase;max_phase]
	obs:select from observations where orbital_phase within (min_phase;max_phase);
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Query by position (circular region on detector)
/ Return all observations where the star centroid falls within a cicular region on the detector
getObsByPosition:{[center_x;center_y;radius]
	/ update adds computed column dist to a new table - doesn't modify global centroids
	/ xexp = power operator, sqrt = square root - implements 2D Pythagorean distance formula across all rows vectorised
	cent:update dist:sqrt[((x_centroid-center_x)xexp 2) + (y_centroid-center_y)xexp 2] from centroids;
	cent_filtered:select from cent where dist<=radius;				/ Filter to rows within circular radius
	0!observations lj `obs_id xkey cent_filtered lj `obs_id xkey psf_params		/ Join filtered centorids and psf onto full observations table
	};

// Get all data for a session
/ Return all observations belonging to s specific session (exact foreign key match)
getSessionData:{[sess_id]
	obs:select from observations where session_id=sess_id;		/ = exact match rather than within range
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Find transit events (flux drops below threshold)
/ Return observations in time range where flux drops below threshol - identifies planetary transit events
getTransitEvents:{[start_bjd;end_bjd;max_flux]
	obs:select from observations where BJD_TDB within (start_bjd;end_bjd), flux<max_flux;	/ comma = AND -> both conditions must be true
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Find high amplitude events in time range
/ Return high brightness amplitude events in a time range
getHighAmpEvents:{[start_bjd;end_bjd;min_amp]
	obs:select from observations where BJD_TDB within (start_bjd;end_bjd), amplitude>min_amp;
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
	};

// Find high airmass observations (lower quality data)
/ Return observations with high airmass indicating lower quality data (longer atmospheric path)
getHighAirmass:{[min_airmass]
	obs:select from observations where airmass>min_airmass;		/ Single threshold filter, no range needed
	0!obs lj `obs_id xkey centroids lj `obs_id xkey psf_params
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
	-1 "  Load Log Rows: ",string count target_load_log;
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
-1 "  1. n:loadTargetCSV[\"HATP-32-b.csv\"]";			/ Inferred target from CSV file name
-1 "  2. n:loadTargetCSVAs[\"my_output.csv\";\"HATP-32-B\"]";	/ Explicit target for generic file names
-1 "  3. n:sortTargetInput[\"my_output.csv\";\"HATP-32-B\"]";
-1 "  4. show 5#getObsByTarget[\"HATP-32-B\"]";

-1 "\nQuery Functions:";
/ Documentation strings - / insidde is printed literally, not a comment 
-1 "  sortTargetInput[path; target_name]         / CSV direct; FITS -> pipeline guidance";
-1 "  getObsByTimeRange[start_bjd; end_bjd]     / Query by BJD_TDB";
-1 "  getObsByTarget[target_name]               / Query by target symbol";
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
