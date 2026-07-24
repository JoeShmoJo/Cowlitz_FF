# COWLITZ_INDEX — paste or fetch this first in a new chat

Compact context for the Cowlitz flow frequency work in repo
`JoeShmoJo/Claude`. Read this INSTEAD of crawling the tree; fetch
individual files only as needed. Last updated 24 Jul 2026.

## The three project folders

| Folder | Purpose | Status |
|---|---|---|
| `CAS_Unreg_FF/` | Unreg flow frequency at Castle Rock (USGS 14243000) via MOS holdout — no ResSim | ACTIVE |
| `CAS_Reg_Unreg/` | ResSim input development: inflow back-calc/shaping, locals, volume correction, ensembles | ACTIVE, in development, no memo |
| `Cowlitz_FF_DataPrep/` | 2026 memo-01 / WY_MAX / SSP record | ARCHIVED (also zipped locally) |

## CAS_Unreg_FF method (current)

Hand-clean hourly MOS ELEV iteratively in DSSVue -> official 2014
rating -> hourly storage -> difference = holdout (cfs) -> process
(3-hr roll avg, Oct–Mar, blank sparse-STOR days) -> SSARR-route
(Mayfield + holdout) and (Mayfield alone) to Castle Rock -> routed
difference + observed CAS = hourly unreg. 3/5-day durations: daily
MOS holdout + CAS daily, no routing (lag < 1 day). Gap WYs (bad
hourly ELEV): correct regulated peak to unreg via regression of
(REG − UNREG 1-hr peak) vs max MOS daily ΔS over 1/2/3/4-day windows
near the peak, computed from the clean DAILY elevation record. REG
peak = 1-hr max of USGS hourly (NOT the USGS instantaneous peak
record) so both peaks are computed identically; track reg/unreg
peak-timing offset and exclude pairs > 72 hrs apart.

Run order: `#DataDownload.py` -> `Build_Hourly_Holdout_Unreg.py`
(was Clean_MOS_Holdout.py; smoothing logic removed) ->
`WY_Peak_Records.py` -> `PeakDiff_Storage_Regression.py` (ADOPTED
predictor: dS_2day) -> `Unreg_Durations_MassBalance.py` (3/5-day) ->
`Write_SSP_Record.py` (final assembly -> output/CAS_Unreg_SSP.dss +
wy_record_ssp.csv audit table; pre-reg WY<=1968 uses USGS peak record
+ durations straight from the USGS daily record). QC: `MOS_STOR_RECORD_COUNT.py`,
`MOS_CASTLEROCK_PEAK_DATE_COMPARE.py`. CAS_Unreg_FF runs FIRST; it
owns obsData.dss and the download script.

## Critical facts

- `CAS_Unreg_FF/data/obsData.dss` is the canonical LIVE observed store,
  shared across projects. `//MOS/ELEV//1HOUR/CWMS-CLEAN/` in it is
  HAND-EDITED in DSSVue, irreproducible; some series manually copied
  in. Never regenerate/overwrite casually. CAS_Reg_Unreg scripts read
  and write this file via `../../CAS_Unreg_FF/data/obsData.dss`.
- No record duplication across projects — cross-project reads instead.
- `Cowlitz_FF_DataPrep/data/obsData.dss` (committed) = frozen memo-01
  snapshot; do not confuse with the live store.
- Daily elevation for the regression + 3/5-day mass balance = daily
  mean of hourly CWMS-CLEAN (legacy IRVZZAZD_CLEANED daily record is
  retired; deletion list in CAS_Unreg_FF/docs/DSS_RECORD_CLEANUP.md).
- data/ref_in = source records only; output/ref_out = everything
  scripts write (final products included).
- Planned, not implemented: derive MOS inflows for ResSim from the
  CWMS-CLEAN elevation instead of the raw inflow record.
- Folder structure, both projects: src/, data/, output/,
  diagnostics/, docs/. No ref_data samples; data is committed
  directly. Scripts use a plain REPO_ROOT config (no
  USE_REFERENCE_DATA switch).
- Coding style: flat functions, no classes/argparse, hardcoded paths
  at top, full scripts. DSS: utilsDSS wrapper in
  `CAS_Unreg_FF/src/Cowlitz_Unreg/Cowlitz/` (readDF/writeSeries);
  sentinels <= -900 missing, -902 written for gaps; pydsstools
  conventions per repo STRUCTURE.md history.
- User environment: Windows, network-restricted (SSL inspection),
  runs scripts locally; container/Claude cannot run pydsstools —
  verify logic with stubs/synthetic data instead.

## Open tasks (as of 24 Jul 2026, PM)

1. dS_2day regression ADOPTED. Run Write_SSP_Record.py to produce
   output/CAS_Unreg_SSP.dss; import *-MOS-HOLDOUT pathnames in SSP.
2. Update `CAS_Unreg_FF/docs/MEMO_CAS_Unreg_FF.docx` with the adopted
   dS_2day regression stats and gap-year fills.
3. Prune retired DSS records per CAS_Unreg_FF/docs/DSS_RECORD_CLEANUP.md,
   then squeeze obsData.dss.
4. CAS_Reg_Unreg: implement cleaned-elevation-derived MOS inflows for
   ResSim (future).

## Key file map

- CAS_Unreg_FF: `README.md` (method + run order), `src/*`, output
  `MOS_Cleaned.dss` paths incl.
  `//CASTLE ROCK/FLOW-UNREG//1HOUR/CAS+ROUTED-DIFF/`,
  `wy_peak_records.csv`, `peakdiff_storage_regressions.csv`,
  `unreg_peak_estimates.csv`, `wy_record_ssp.csv`,
  `CAS_Unreg_SSP.dss` (final SSP input), memo in `docs/`.
- CAS_Reg_Unreg: `README.md` (script + diagnostics use cases);
  external ResSim watershed paths under
  `C:\Projects\Cowlitz_Flow_Frequency\ResSim\...` (not in repo).
- Cowlitz_FF_DataPrep: `RUN_ORDER.md` (archive note, diagnostics
  index, memo pipeline), `docs/MEMO_01_...docx`.
