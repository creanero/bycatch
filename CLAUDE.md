# CLAUDE.md

Guidance for Claude Code (and humans) working in this repo.

## What this repo is

Astronomical Bycatch: extract exoplanet transits from stars that happen to be
"in the background" of archival images. See `README.md` for the premise and
`images/structure.svg` for the high-level design (archive → extraction API →
database → short/long-term transit analysis).

Layout:

- `photometry/` — the current photutils-based pipeline (D. Lawley, adapted by
  D. Winters). `target_star_photometry.py` aligns frames with astroalign and
  does target/comparison photometry; `bycatch_photometry.py` detects every
  source with DAOStarFinder and photometers them all; `fwhm_calibration.py`
  derives FWHM and aperture/annulus radii from the reference frame;
  `config.py` holds paths and parameters (still hard-coded, Windows paths).
- `database/` — kdb+/q schema (`targetStarOnly.q`: `sources`,
  `obs_sessions`, `observations`) and its design notes.
- `publications/`, `images/` — posters and figures.

## Relationship to EXOTIC (`~/work/exotic-gitlab/exotic-private`)

This repo is the consumer of the new extraction API being built in the
private EXOTIC fork (GitLab canonical; see that repo's `CLAUDE.md` and
`docs/tests-and-extraction-api.md`). **The intention is to use that API to get
light curves for the targets we define here**, rather than growing a second,
parallel photometry stack in `photometry/`.

The API lives in `exotic/extract/` (not `exotic/api/`, which is EXOTIC's
older library code). It does no file I/O and no source search — the caller
reads frames and supplies positions.

- **Layer 0** primitives: `photometry.py` (`fit_centroid`, `skybg_phot`,
  `aperture_sum`, `grid_photometry`), `timing.py` (JD/BJD_TDB, airmass),
  `calibration.py` (`apply_cals`), `alignment.py` (`find_transform`,
  astroalign-based), `uncertainty.py` (`ccd_flux_uncertainty`).
- **Layer 1** `frame.extract_frame(frame_data, positions, *, aperture_r,
  annulus_r_in, annulus_r_out, gain, read_noise, bias_level_per_pixel=0)` →
  list of `ExtractedSource` (centroid, flux, uncertainty, PSF params). A
  failed centroid gives a NaN row rather than raising.
- **Layer 2** `sequence.extract_sequence(frames, positions=None, *,
  sky_positions=None, ra, dec, observatory, apertures, gain, read_noise,
  bias_level_per_pixel=0, dark=None, bias=None, flat=None, align=True)` →
  numpy structured array (`SEQUENCE_DTYPE`), one row per frame per source,
  frame-major, with `bjd_tdb`, `jd_utc`, `airmass`, `aligned`, the radii, the
  `ExtractedSource` fields and the centroid's `ra`/`dec`.
  - `frames` is an iterable of `(data, header)`; the first is the reference.
    A generator keeps memory flat.
  - Sources as reference-frame pixels or as RA/Dec (`sky_positions`, via the
    reference header's WCS). RA/Dec out lets the same star be linked across
    sequences — this is what the database needs.
  - Bad frames are flagged, not dropped; mask on `aligned` and finite `flux`.
- **Aperture strategies** (`apertures.py`): `Apertures(r, r_in, r_out)`
  fixed, or `FwhmScaled()` which measures FWHM on the reference frame. Its
  defaults (1.9, 1.1, 2.1) reproduce this repo's `config.py` coefficients
  (2.6, 1.55, 2.9) on HAT-P-32; the difference is only that EXOTIC's Gaussian
  FWHM is ~1.38x IRAFStarFinder's. On new data they are a starting point.
- **Transit fit**: `exotic.api.quick_fit.fit_transit(...)`.
- Worked examples: `examples/extract_api/{single_frame_extract,
  light_curve_extract, light_curve_fit}.py` in the EXOTIC repo.

Known open items on the EXOTIC side that affect us:

- Gain / read noise / bias pedestal are not yet calibrated for any dataset;
  the examples use placeholders, so uncertainties are not yet scientific.
  The HAT-P-32 validation is blocked on getting real dark frames.
- astroalign's RANSAC is unseeded, so alignment (and hence centroids) is not
  bit-reproducible run to run.

EXOTIC runs inside its Docker image (`make` targets in that repo), not host
Python; use that environment when calling the API from here.

## Target data

Still being defined. Targets we pick should be recorded here (or in a data
file this section points to) with: field RA/Dec, source list (RA/Dec),
archive and query used, instrument, observatory site, and which calibration
frames exist.

### Data source: ESO Science Archive

<https://archive.eso.org/eso/eso_archive_main.html> — the ESO raw-data query
form for La Silla Paranal instruments. Notes so far:

- Searches **raw frames** only (reduced Phase 3 products are a separate
  portal). Queries by target (SIMBAD/NED resolved), RA/Dec box, night or
  date range, programme ID/type, instrument, exposure time, filter, and by
  category (SCIENCE/CALIB/…), type (OBJECT/BIAS/DARK/FLAT…) and mode (IMAGE,
  SPECTRUM, …). Output as HTML, ASCII or VOTable; downloads can include
  associated calibrations.
- Programmatic access (verified): the ESO TAP service handles metadata
  search; `astroquery.eso.Eso().retrieve_data([dp_id])` downloads public
  frames without a login (it fetches `.fits.Z` and decompresses it).
- Why it is attractive for bycatch: long imaging time series from many
  programmes, with **bias/dark/flat calibrations archived alongside** — which
  is exactly what EXOTIC's flux-uncertainty work is currently missing.
- Things to check before committing to it:
  - Proprietary periods: only data past its proprietary period is public.
  - Which imaging instruments/modes give useful bycatch (field of view, depth,
    cadence, sequence length) — restrict to `mode = IMAGE`, science `OBJECT`
    frames, and long same-field sequences in one night.
  - Raw ESO frames are often multi-extension FITS (mosaic detectors); each
    detector/extension is its own `(data, header)` sequence for
    `extract_sequence`.
  - Header keywords: ESO uses `HIERARCH ESO ...` keywords; confirm EXOTIC's
    `timing.py` reads the right time/exposure keywords, and that the raw WCS
    is good enough for `sky_positions` (or plate-solve first).
  - Observatory site (`Observatory(lat, lon, elev)`) per telescope —
    Paranal vs La Silla — ideally read from the header's telescope-site
    keywords rather than hard-coded.
  - Gain / read noise per instrument and readout mode, from headers or
    instrument documentation.

### ESO survey results (20 Oct 2026, metadata only)

Queried via the ESO TAP service (`https://archive.eso.org/tap_obs`,
table `dbo.raw`) and the NASA Exoplanet Archive TAP (`pscomppars`,
`tran_flag=1`). Scripts were ad hoc; turn them into `access/` code before
relying on them.

- Public science imaging frames: OmegaCAM ~447k, VIRCAM ~2.0M (NIR,
  jittered pawprints, poor for transits), WFI ~189k, FORS2 ~264k, HAWK-I
  ~189k, SOFI ~792k, EFOSC ~161k.
- ESO's ADQL rejects `FLOOR` / grouping by expressions, so nights have to be
  formed client side: pull per-frame metadata (`instrument, prog_id, object,
  filter_path, mjd_obs, ra, dec, exposure`) and group by
  `night = floor(mjd_obs - 0.5)`. Pulling all WFI + OmegaCAM rows (~640k)
  takes a few minutes with `run_async`.
- Long same-pointing sequences (≥40 frames, ≥2 h, pointing scatter <0.1°):
  132 WFI, 50 OmegaCAM. Only two known transiting planets fall in any of
  them; only one has a transit in the window.

**Test target: POTS-1 b, ESO/MPG 2.2 m WFI, night of 2007-05-10/11.**

- Field `OTSF-1A` (Pre-OmegaTranS), programme `079.A-9003(A)`, Rc filter
  (`BB#RC/162_ESO844`), 310 frames × 26 s over 8.2 h. Ephemeris predicts a full
  transit with 2.5 h of baseline before and 2.8 h after (T_c uncertainty
  <1 min). POTS-1 is V≈17.9 with a deep transit; this survey is where it was
  found, so it tests whether we recover a known result, not an independent
  detection. The field is dense Galactic disk, so there's plenty of bycatch.
- Calibrations the same night: 10 Rc biases, 14 Rc sky flats, 5 Rc screen
  flats.
- Size: ~60 MB compressed / 142 MB uncompressed per frame, so ~18 GB
  compressed for the night.
- One frame (`WFI.2007-05-11T03:46:56.142`) inspected. What it shows:
  - 8 image extensions (`WIN1.CHIPn.OUTm`), each 4128×2142 *including
    overscan/prescan*, which must be trimmed before photometry.
  - **No celestial WCS** (`CTYPE = PIXEL`). `sky_positions` won't work on
    raw frames: plate-solve the reference frame (per chip) first.
  - Site in `HIERARCH ESO TEL GEOLAT/GEOLON/GEOELEV` (−29.2543, −70.7346,
    2335 m); time in `MJD-OBS`/`DATE-OBS`, `EXPTIME`; airmass in
    `ESO TEL AIRM START/END` (there's no plain `AIRMASS` keyword).
  - Header gain/read noise are unreliable (0.0 or missing on most chips), so
    take them from the WFI instrument documentation.

### Data source: MicroObservatory frames (collaborator's private repo)

Raw MicroObservatory frames of several hot-Jupiter targets (WASP-36,
HAT-P-32, HAT-P-36, WASP-12, Qatar-2, HAT-P-54, WASP-140). These are
target-centred frames from the camera EXOTIC is already validated on, so
they suit per-target extraction checks and bycatch of the field stars. ESO
WFI is the wide-field, deep complement. Where the repo is, and where it is
cloned: `CLAUDE.local.md`.

## Plan

1. **Data access (`access/`)**: astroquery/pyvo search → sequence table →
   download. Keep data on the compute host's data share, not on
   laptops. Stream it: a generator downloads a frame, yields one chip's
   trimmed `(data, header)` to `extract_sequence`, and the cache can be
   cleared afterwards. Plan for one full download (~18 GB) for the test
   night, since `extract_sequence` needs one pass per chip.
2. **Instrument adapter**: per-instrument trimming, chip layout,
   gain/read noise, site and time keywords, and building master bias/flat
   from the same night's CALIB frames.
3. **Astrometry**: no local astrometry.net on the compute host. Options are a Gaia
   query (`astroquery.gaia`, available) matched to detected sources, giving
   a fitted WCS per chip of the reference frame, or installing
   astrometry.net with the index files for a ~0.5° field. Prefer Gaia
   matching: no external upload, and Gaia IDs give the database a stable
   source key.
4. **Extraction**: `exotic.extract.extract_sequence` with Gaia
   `sky_positions` (target + all detected sources), `FwhmScaled()`, real
   calibrations. Run on the compute host.
5. **Validation**: `exotic.api.quick_fit.fit_transit` on POTS-1, comparing
   T_mid and depth with the published ephemeris.
6. **Results store (JSON for now)**: one file per (sequence, chip), e.g.
   `results/<instrument>/<night>/<chip>.json`, holding sequence metadata
   (archive, dp_ids, instrument, filter, site, apertures, gain, exotic
   version) and per-source light curves keyed by Gaia DR3 ID (ra, dec,
   bjd_tdb[], flux[], uncertainty[], flags). Maps onto the q schema's
   `sources` / `obs_sessions` / `observations` later.

## Data access code (`access/`, `targets/`)

- `access/eso.py`: TAP queries (`query_frames`, `find_sequences`,
  `sequence_frames`, `night_calibrations`) and a threaded `download` that
  skips files already present.
- `access/transits.py`: NASA Exoplanet Archive ephemerides (`transiting_planets`)
  and `transits_in(sequences, planets)` to find sequences catching a transit.
- `access/fetch.py`: `python -m access.fetch targets/<t>.json [--data-root
  ~/bycatch-data] [--manifest-only]`. Writes
  `<root>/eso/<instrument>/<night>_<object>/manifest.json` (every science and
  calibration frame's ESO metadata) and downloads into `science/` and
  `calib/`. Restartable.
- `targets/*.json`: one per target sequence. Gives archive, instrument,
  prog_id, object, filter, night (`floor(MJD-0.5)`) and `calibrations`:
  dp_type → `true` (same filter as the science frames) / `false` (any filter).
- Data root defaults to `~/bycatch-data`. The POTS-1 test night is
  `<root>/eso/wfi/54230_OTSF-1A/`: 310 science + 38 calibration frames,
  ~50 GB uncompressed.

## Pipeline code (`pipeline/`)

- `pipeline/astrometry.py`: local plate solving against Gaia DR3 (DAOStarFinder
  detections + astroalign pattern match + `fit_wcs_from_points` with SIP).
  `solve_wcs_robust` tries several centres and settings. Needed because raw
  WFI and MicroObservatory frames have no WCS, and MicroObservatory header
  pointing can be ~0.3 deg off.
- `pipeline/store.py`: JSON results store, schema `bycatch-lightcurves/1`
  (sequence metadata, per-frame times/airmass, per-source light curves
  keyed by Gaia DR3 ID, analysis).
- `pipeline/microobs.py`: MicroObservatory night → solve → Gaia sources →
  `extract_sequence` → comparison ensemble → `fit_transit` with airmass
  detrending → JSON + PNG. `python -m pipeline.microobs FRAMES_DIR --inits
  INITS.json --out OUT.json [--plot OUT.png] [--aperture-coeff C]`. Needs
  EXOTIC on `PYTHONPATH`.

First run, WASP-36 b (MicroObservatory, 2026-01-05, 55 × 60 s):

- Solve: 91 Gaia stars, 0.94" RMS, 5.237"/px. Tracking holds: target vs.
  comparison offset scatter <0.6 px across ~70 px of drift.
- EXOTIC's JD from the header equals the DATE-OBS-derived mid-exposure JD
  exactly on these frames, so there's no MJD-OBS rounding problem here.
- Noise-limited: per-point scatter goes from 1.6% at G≈11 to 6% at the target
  (G=12.6) with the default aperture (1.9 × FWHM). A sweep over
  `--aperture-coeff` 1.0–1.9 gives monotonically lower scatter for smaller
  apertures: 2.3% at 1.0 × FWHM (2.3 px). That suggests sky-dominated noise,
  and is worth checking against the background estimate.
- Transit (~1.8% deep) only weakly detected. Best fit (1.0 × FWHM): O−C
  +26 ± 43 min, Rp/Rs 0.09 ± 0.06, vs. predicted T_mid 2461045.9342 ±
  0.4 min. A 66 min gap leaves only 3 frames of pre-ingress baseline, and
  there are no flats or darks.

Validation run, HAT-P-32 b (MicroObservatory, 2017-12-20, 142 × 60 s; the
same frames as EXOTIC's `HatP32Dec202017` sample data):

- Solve: 184 Gaia stars, 1.5" RMS. Target G=11.1; 83 sources photometered.
- Against EXOTIC's full-reduction validation anchor (2458107.7137986
  BJD_TDB):

  | aperture | radius | scatter | T_mid − anchor | Rp/Rs |
  | --- | --- | --- | --- | --- |
  | 1.0 × FWHM | 2.5 px | 0.53% | −0.6 ± 2.6 min | 0.151 ± 0.003 |
  | 1.4 × FWHM | 3.5 px | 0.48% | +2.1 ± 2.4 min | 0.157 ± 0.003 |
  | 1.9 × FWHM | 4.7 px | 0.57% | +4.0 ± 2.8 min | 0.156 ± 0.004 |

  All within ~1.5σ of the anchor. Published Rp/Rs is 0.149. The fit is a
  least-squares `fit_transit` with fixed rough limb darkening and a linear
  airmass trend, not EXOTIC's nested-sampling fit, so a few minutes'
  aperture-dependent spread is expected.
- The optimum aperture is ~1.4 × FWHM here, not the smallest. WASP-36's
  monotonic trend is a faint-star (sky-limited) effect, not a pipeline bug.
- Inits coordinates can be rounded by ~10"; the target is matched to the
  brightest Gaia star within 15".

WFI POTS-1 b (night 54230, chip 3, 801 sources, `python -m pipeline.wfi
calibrate|solve|extract`): recovered at 1.2 × FWHM with O−C +4.7 ± 4.0 min,
Rp/Rs 0.166 ± 0.010 (catalogue 0.164), 2.3% rms. The report is
`reports/progress_2026-10-20.tex`.

POTS-1 b second night (54173, `targets/pots1_wfi_54173.json`): O−C +1.3 ±
2.2 min, Rp/Rs 0.165 ± 0.007, 1.3% rms. With 54230 (18 orbits apart):
P = 3.16076 ± 0.00018 d vs catalogue 3.16063 d (0.7σ).

WASP-4 b, EFOSC2 (`pipeline/efosc.py`, `targets/wasp4_efosc_55855.json`):
O−C +0.19 ± 0.29 min against a ±0.02 min prediction; Rp/Rs 0.1518 ±
0.0007 (catalogue 0.1520); rms 0.12% at 1.9 × FWHM. EFOSC2 specifics:

- Windowed readout (700 × 900 at 0.12"/px), with full-frame calibrations.
- Header WCS is in full-frame pixels.
- Defocused PSF (~30 px), so FWHM comes from the radial profile.
- Only 3–4 stars, so astroalign fails. Registration tracks the target with
  whole-pixel shifts.
- Sky flats are `SKY,FLAT`, not WFI's `FLAT,SKY`.
- Acquisition frames mixed in, so the target file uses `"exposure"`.

Wider ESO search (7 instruments, cached under `~/bycatch-data/search/` on
cheetah): 92 sequence/transit matches for 41 planets. Mostly target-centred;
OGLE-TR fields (FORS1/2, SOFI) are the crowded, bycatch-like ones.

Archive outages seen: the NASA Exoplanet Archive returned HTTP 526 for over
an hour, and Gaia gave intermittent 500s. Planet rows are cached in
`~/bycatch-data/cache/planets/`; Gaia queries retry.

Blind searches (no knowledge of the target):

- `pipeline/search.py` (one night). On POTS-1 chip 3, night 54230, it
  ranks POTS-1 9th of 764 (99th percentile), the only transit-shaped event
  in the top 12.
- `pipeline/periodsearch.py` (several nights, stars matched by Gaia ID,
  BLS). With nights 54173 + 54230, POTS-1 is #1 of 486 multi-night
  solutions at SNR 23.1; the next is 6.7. It passes strict vetting:
  - each event significant on its own (2.87 ± 0.11% and 2.28 ± 0.22%);
  - depths consistent;
  - both events complete within their nights.

  The period is ambiguous with two nights. The catalogue 3.16066 d is
  among 76 aliases of the 56.894 d gap.

**Whole POTS-1 field, both nights** (22 Oct 2026, overnight on the compute host):

- Run with `scripts/run_pots_field.sh` on the compute host. Field-only `pipeline.wfi
  extract --chip N` for the 7 other chips, about 11 min per chip on the
  short night and 21 min on the long one, with 16 workers. Then
  `periodsearch` per chip and `pipeline.collate`. Results are under
  `~/bycatch-data/results/eso/wfi/pots_field/` (`field_candidates.json`,
  `field_top8.png`).
- About 5,900 stars searched on both nights.
- POTS-1 is #1 field-wide (SNR 23.1; next 18.8). By eye it's the only
  convincing transit. The next candidates are noise, trends or one-sided
  events. The only one worth a second look is #5 (chip 2, G=14.8, Gaia
  5851641466227710208), a weak 0.7% box on both nights.
- 126 stars excluded as variable. 33 have a Gaia DR3 class (19 ECL, 6
  DSCT/GDOR, 5 RS, 1 RR, among others). The other 93 are unclassified:
  possible new variable-star bycatch.
- Weaknesses to fix:
  - "Complete" events should need ≥30 min of data on both sides (the #4
    candidate's event starts at the night's first frame).
  - Two Gaia LPVs slipped past the variable flag.
  - Chip 8 is noisier than the rest, so noise should be judged per chip.
  - Weak candidates pile up near the 0.7 d minimum period.

What made the difference, in order:

1. Cotrending on the field's principal components (not a single common
   mode).
2. Clipping single-frame outliers.
3. Scoring against noise on the transit timescale, estimated from sliding
   windows.
4. Excluding variables by the spread left after removing the best box (a
   deep transit on a short night otherwise looks variable).
5. A field-median correlated-noise factor (a star's own factor is inflated
   by its transit).
6. Masking periods near n days and 1/n day (time-of-night systematics).
7. Vetting per event.

Remaining bias: weak candidates pile up near the 0.7 d minimum period.

EXOTIC pitfalls found here, now fixed in exotic-private `main` (801067c).
The bycatch workarounds are removed. The batched grid fix is also on `main`
(0cba1cc), and the compute host runs a git clone of `main` (see
`CLAUDE.local.md`). Make sure any exotic copy includes 801067c.

- `skybg_phot` used to take the *mode* of a box around each source, which
  on flat-fielded float data is the smallest pixel (20–30% flux errors).
  On non-integer data it now uses a clipped annulus median, so calibrated
  frames go in as floats. That took WASP-4 from 0.12% to 0.05% rms and
  POTS-1 from 2.3% to 1.1%.
- `quick_fit.fit_transit` used to fit absolute BJD, where curve_fit's
  finite-difference step on tmid is ~50 min, so fits stalled. It now fits
  relative to the night. `lightcurve.fit_with_airmass` still adds a grid
  start and 4σ clipping. `python -m pipeline.refit RESULTS.json` refits
  stored results.
- Centroid fitting costs ~70 ms per source per frame. `pipeline.wfi`
  splits sources across processes (`--workers`).

## Local setup

Machine-specific details (compute host, environments, data locations,
private data sources) are in `CLAUDE.local.md`. It is gitignored: keep
hostnames, accounts and private repo details there, not here, because this
repo is public.

## Working conventions

- Branch names include the author's name, e.g. `huw-leggate/<topic>`.
- Don't commit bytecode (`__pycache__`) or large data/FITS files.
