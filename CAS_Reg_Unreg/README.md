# CAS_Reg_Unreg

ResSim input development for the Cowlitz regulated/unregulated work:
reservoir inflow back-calculation and shaping, local flow calculation,
volume correction, and ensemble DSS creation. Process still in
development -- no memo yet.

Kept separate from `CAS_Unreg_FF/` (the Castle Rock unreg flow
frequency record, which no longer depends on ResSim). Scripts here may
READ shared records from `CAS_Unreg_FF/data/obsData.dss` (the canonical
live observed-data store) rather than duplicating them; several also
WRITE derived records (calculated inflows, volume-corrected series,
locals) back into that file.

## Scripts (src/)

Data download lives in CAS_Unreg_FF (`#DataDownload.py`) -- that
project runs first and owns the shared obsData.dss.

- `#Reservoir_Inflow_And_Basin_Peaks.py` -- inflow back-calc from
  elev/stor + outflow, peak-window shaping, Castle Rock local flow
  inputs. WRITES into shared obsData.dss; do not run casually.
- `#Inflow_Volume_Correction.py` -- scales cleaned/smoothed inflow
  segments to match raw calculated volumes (FLOW-IN-CALC-CLEANED-
  VOLCOR); also shapes the Mayfield local.
- `#Create_Ensembles.py` -- builds ensemble.dss for the ResSim
  watershed (external path) + mapping CSV; computes Castle Rock local
  (CAS - Mayfield - Toutle).
- `#ExtractResSimEnsembleResults.py` -- pulls ResSim simulation output
  (external path) -> `data/results.dss` (WY 1/24/72/120-hr maxima).
- `#MOS_CDB_INFLOW.py` -- 1996 MOS inflow from CDB elevation via mass
  balance (ESRD testing input).

## Diagnostics (diagnostics/) -- use cases

- `*_peak_windows.html` -- inspect shaped inflow around event peaks
  (MAY, MOS) before accepting the shaping.
- `*_inflow_comparison.*` -- raw vs cleaned vs volume-corrected inflow
  overlay per reservoir.
- `MOS_volcor_comparison.html`, `MOS_volcor_stats.csv`,
  `inspect_VOLCOR.csv` -- volume-correction segment scales and QA.
- `MAY_local_shaped.html`, `MAY_local_shaped_stats.csv`,
  `inspect_MAY_PEAKCLEAN.csv` -- Mayfield local shaping QA.
- `2015/2025_MOS_INFLOW_CLEAN.csv` -- example cleaned-inflow exports.
- `alignment_check.csv`, `ensemble_slot0_check.csv` -- ensemble
  alignment verification.

## Planned change (not yet implemented)

MOS inflows for ResSim should eventually be calculated from the
hand-cleaned elevation record in CAS_Unreg_FF (CWMS-CLEAN) instead of
the raw inflow record used previously.
