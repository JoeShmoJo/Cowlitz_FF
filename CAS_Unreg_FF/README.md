# CAS_Unreg_FF

Unregulated flow frequency record at Castle Rock (USGS 14243000) built
by the MOS holdout method -- NOT by routing cleaned inflows through
ResSim. Formerly `Basic_CAS_Unreg`; ResSim input development lives in
`CAS_Reg_Unreg/`, and the 2026 memo-01 record work is archived in
`Cowlitz_FF_DataPrep/`.

## Method in one paragraph

The hourly MOS elevation record is hand-cleaned iteratively in DSSVue
(edit -> compute holdout -> scrutinize -> edit again). The cleaned
elevation is converted to storage with the official 2014 rating and
differenced to a raw hourly holdout, which is processed (3-hr rolling
average, Oct-Mar season, sparse-telemetry days blanked) and routed with
Mayfield outflow to Castle Rock via three chained SSARR reaches. The
routed holdout effect added to observed Castle Rock flow gives the
hourly unregulated estimate. For 3- and 5-day durations, daily holdout
at MOS added to Castle Rock daily flow is adequate (routing < 1 day,
so no routing needed). For WYs where the hourly elevation was too poor
for a holdout, the peak regulated flow is corrected to unregulated via
a regression of (REG - UNREG peak) against daily MOS storage change.

## Manual step -- do not lose this record

`//MOS/ELEV//1HOUR/CWMS-CLEAN/` in `data/obsData.dss` is hand-edited
in DSSVue and is NOT reproducible from any script. It is the anchor of
the whole workflow. Some series were also manually copied from the
observed data into this file. Backup accordingly.

## Scripts (src/), in run order

0. `#DataDownload.py` -- USGS/CWMS download -> `data/obsData.dss` (+
   optional parquet). Record list read from `data/MOS_ELEV.csv`
   (RequiredRecordsDict format; fuller dictionaries alongside).
   CAUTION: never let a re-download overwrite the hand-cleaned
   CWMS-CLEAN record -- it writes the raw CWMS/USGS records only, but
   back up obsData.dss first anyway.
1. `Build_Hourly_Holdout_Unreg.py` (formerly Clean_MOS_Holdout.py --
   renamed because the despike/smoothing machinery was removed once
   ELEV moved to hand-cleaning). Cleaned ELEV -> holdout -> routing ->
   hourly unreg. Writes `MOS_Cleaned.dss`.
2. `WY_Peak_Records.py` -- per-WY 1-hr and 1-day peaks of the hourly
   unreg and hourly regulated (1-hr max of USGS hourly, NOT the USGS
   instantaneous peak record), peak timestamps, REG-UNREG diffs, and
   the timing offset between the two peaks. -> `wy_peak_records.csv`.
   Also records unreg_cov_at_reg_peak (unregulated data availability
   within +/-24 hrs of the regulated annual peak) -- the screen the
   assembly uses to accept or reject an hourly peak -- and identifies
   every Oct-Mar missing window >= MIN_GAP_HRS in both series (-> `diagnostics/wy_missing_windows.csv`, with each gap's
   distance to that WY's peak, plus per-WY gap-summary columns in the
   main table). Nothing is omitted automatically -- review the gap
   report, then populate EXCLUDE_RANGES (mask a suspect window) or
   OMIT_WYS (drop a whole WY) and re-run.
3. `PeakDiff_Storage_Regression.py` -- regress (REG - UNREG) 1-hr peak
   against max MOS daily storage change over 1/2/3/4-day windows near
   the peak (daily means of CWMS-CLEAN). ADOPTED PREDICTOR: dS_2day.
   Applies the fit to gap WYs that have a good regulated peak (hourly
   or USGS peak record) but no holdout-based unreg peak.
4. `Unreg_Durations_MassBalance.py` -- daily mass-balance unreg
   (CAS daily + dS from the daily MOS elevation record
   //MOS/ELEV//1DAY/USGS/), 1/3/5-day WY maxima. Used for the 3- and
   5-day durations; no routing needed at daily resolution. The daily
   record is separate from -- not derived from -- the hourly
   CWMS-CLEAN, which is used only in the holdout workflow under its
   STOR-COUNT day-by-day screen. Admission: Season_Complete =
   Oct-Mar completeness screen OR SEASON_OVERRIDE_WYS (adopted:
   the full regulated era WY1974-2025 -- the Castle Rock daily record
   is winter-only / season-incomplete through much of it -- under the
   stated assumption that the recorded season captured the annual
   maxima; WY2026 excluded until the year closes). Rolling windows run
   on a full daily grid, so an N-day value always means N consecutive
   recorded days.
5. `Write_SSP_Record.py` -- assemble the final Peak/1/3/5-day WY record
   and write `output/CAS_Unreg_SSP.dss` for HEC-SSP, plus the audit
   table `wy_record_ssp.csv` (every WY carries Peak_Screen, both peak
   candidates, and the coverage evidence behind the choice) and the QA
   file `diagnostics/record_qa_flags.csv` (duration-ordering
   violations, classed SAME_SOURCE vs cross_source, logged BEFORE any
   adjustment). With ENFORCE_MONOTONIC (default on) the record is then
   made monotonic working backwards from 5-day: Five_Day is the anchor,
   each shorter duration is raised to the longer one where the longer
   is higher, and raises cascade upward. Pre-adjustment values are kept
   in the *_Raw columns, each WY's changes are summarized in
   Monotonic_Adjustment, and every raise is logged to
   `diagnostics/record_monotonic_adjustments.csv`. Missing shorter
   durations are left missing, never fabricated from a longer one.
   EXCLUDE_WYS drops a whole water year from the record whatever source
   would have supplied it, with the reason and the discarded values
   written to `diagnostics/record_excluded_wys.csv`. (Note:
   SEASON_OVERRIDE_WYS in the mass-balance script governs ONLY the
   daily durations -- removing a WY there still leaves its Peak and
   hourly One_day in the record.) Source rules: pre-1968, Peak from the
   USGS peak flow record and 1/3/5-day straight from the USGS daily
   record. Post-1968, Peak from the calculated hourly unreg else the
   dS_2day storage-change regression; One_day from the 1-day average
   of the hourly unreg else the one-day max of the daily unreg; 3- and
   5-day from the unreg daily averages (mass balance). Every value is
   source-tagged. Nothing reads from the Cowlitz_FF_DataPrep archive
   (CastleRock_USGS_peaks.csv now lives in data/).

QC / reference:
- `MOS_STOR_RECORD_COUNT.py` -- daily count of valid hourly STOR values
  (spot sparse-telemetry days before trusting a holdout).
- `Cowlitz_Unreg/Cowlitz/` -- utilsDSS wrapper (readDF/writeSeries,
  sentinel handling), SSARR routing, config.

## Folder structure (consistent across projects)

`src/`, `data/`, `output/`, `diagnostics/`, `docs/`. No ref_data
sample folders -- the data is small enough to commit directly.
Input/output rule: `data/` holds SOURCE records only (observed
downloads + the hand-cleaned store); `output/` holds everything
scripts write, intermediates (MOS_Cleaned.dss) and final products
(CAS_Unreg_SSP.dss) alike. See docs/DSS_RECORD_CLEANUP.md for the
record-level cleanup list and the one Reg_Unreg exception.
- `data/obsData.dss` is the canonical live observed-data store, shared
  with CAS_Reg_Unreg (its scripts point here). Don't duplicate records
  across projects -- cross-project reads are preferred.
- Sentinels <= -900 are missing; -902 written for gaps.

## Status (24 Jul 2026): RECORD COMPLETE

WY1927-2026 assembled and written to output/CAS_Unreg_SSP.dss; only
WY1969-1973 absent (perception-threshold candidates in SSP). Peak: 95
values (42 USGS pre-reg / 37 dS_2day regression / 16 hourly holdout);
One_day 94; Three/Five_Day 93 (WY1927 durations n/a; WY2026 3/5-day
await year completion). Adopted regression: (REG - UNREG) =
-0.869 * dS_2day - 3,836 cfs, R^2 0.871, SE ~5,530 cfs, n 17. Memo
finalized in docs/. Remaining: WY1969-1973 SSP treatment; low-coverage
fit-set review; WY2026 close-out re-run; DSS cleanup + squeeze per
docs/DSS_RECORD_CLEANUP.md.
