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
   Also identifies every Oct-Mar missing window >= MIN_GAP_HRS in both
   series (-> `diagnostics/wy_missing_windows.csv`, with each gap's
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
   (CAS + daily MOS holdout), 1/3/5-day WY maxima. Used for the 3- and
   5-day durations; no routing needed at daily resolution.
5. `Write_SSP_Record.py` -- assemble the final Peak/1/3/5-day WY record
   and write `output/CAS_Unreg_SSP.dss` for HEC-SSP, plus the audit
   table `wy_record_ssp.csv`. Pre-regulation WYs (<= 1968, before
   Mossyrock closure): peaks from the USGS instantaneous peak record,
   1/3/5-day computed directly from the USGS daily record (both
   unregulated by definition pre-dam). Regulated era: hourly peaks +
   dS_2day regression fills + mass-balance durations. Every value is
   source-tagged.

QC / reference:
- `MOS_STOR_RECORD_COUNT.py` -- daily count of valid hourly STOR values
  (spot sparse-telemetry days before trusting a holdout).
- `MOS_CASTLEROCK_PEAK_DATE_COMPARE.py` -- peak-date comparison
  diagnostics between MOS storage and Castle Rock records (moved from
  Cowlitz_FF_DataPrep; reads that archive's data).
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

## Status / open items

- Hourly unreg constructed; iterate ELEV cleaning as needed.
- Run WY_Peak_Records + PeakDiff_Storage_Regression on real data;
  evaluate which dS window wins; scrutinize WYs with large peak-timing
  offsets before accepting them into the fit.
- Brief memo in docs/ describes the methodology; update with regression
  results once adopted.
- Prune retired DSS records per docs/DSS_RECORD_CLEANUP.md, then
  squeeze obsData.dss.
