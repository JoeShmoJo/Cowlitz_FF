# DSS Record Cleanup — obsData.dss (live store)

Compiled 24 Jul 2026 by cross-referencing every DSS pathname referenced
in the active scripts (CAS_Unreg_FF/src + CAS_Reg_Unreg/src) against
the record inventory of the obsData.dss sample. Delete in DSSVue on the
LIVE file only. D parts omitted (all blocks of each family). The frozen
copy in Cowlitz_FF_DataPrep/data/obsData.dss and the local zip archive
retain everything, so nothing below is lost to history.

After deleting, run a Squeeze in DSSVue — DSS6 keeps deleted-record
space (and remnants of old renames) until the file is squeezed.

## KEEP — actively used

| Record | Used by |
|---|---|
| //MOS/ELEV//1HOUR/CWMS-CLEAN/ | Build_Hourly_Holdout_Unreg (gated day-by-day by the STOR-COUNT screen + overrides). HAND-EDITED, irreplaceable. |
| //MOS/ELEV//1DAY/USGS/ | Daily MOS elevation -- Unreg_Durations_MassBalance and PeakDiff_Storage_Regression dS predictors |
| //MOS/STOR//1HOUR/CWMS/ | Holdout QC count, MOS_STOR_RECORD_COUNT, Reg_Unreg inflow calc |
| //MOS/ELEV//1HOUR/CWMS/ | Raw source the hand-cleaning started from; re-downloadable but cheap provenance — keep |
| /COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW//1HOUR/USGS/ | Holdout, WY_Peak_Records, Reg_Unreg |
| /COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW//1DAY/USGS/ | Unreg_Durations_MassBalance |
| /COWLITZ RIVER BELOW MAYFIELD DAM, WA/14238000/FLOW//1HOUR/USGS/ | Holdout routing, Reg_Unreg |
| /COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW-LOCAL//1HOUR/USGS/ | Reg_Unreg |
| Packwood 14226500, Cispus 14232500, Tilton 14236200, Toutle 14242580 FLOW records | Reg_Unreg (download dict, ensembles) |
| //MOS + //MAY FLOW-IN / FLOW-OUT / STOR //1HOUR/CWMS/ families | Reg_Unreg inflow calc |
| //MOS + //MAY FLOW-IN-CALC-RAW / -RAW-PEAKS / -CLEANED / -CLEANED-VOLCOR | Reg_Unreg (derived, regenerable, but feed ResSim) |
| //MAY/FLOW-LOCAL, FLOW-LOCAL-SHAPED, FLOW-OUT_PEAKCLEAN_2009_2026, STOR_PEAKCLEAN_2009_2026 | Reg_Unreg |
| //MOS/ELEV-RULECURVE//1HOUR/CENWP-CALC/ | #ExtractResSimEnsembleResults |
| /ZERO_DUMMY/ZERO/FLOW//1HOUR/DUMMY/ | #Create_Ensembles |

## DELETE — no active script references them

| Record family | Why it existed / why it can go |
|---|---|
| //MOS/ELEV-FOREBAY//1DAY/IRVZZAZD_CLEANED/ (+ //MAY/... if present) | Legacy daily cleaned elevation. Confirmed retired; scripts now derive daily means from CWMS-CLEAN. Reproducing archived memo results uses the frozen dataprep obsData.dss. |
| //MOS/FLOW-HOLDOUT//1HOUR/RAW/, /RAW_3H_AVE/, //MOS/FLOW-HOLDOUT//IR-YEAR/RAW/ | Early holdout iterations written into the input store. Current outputs go to MOS_Cleaned.dss. |
| //MOS/FLOW//1HOUR/CWMS_WINTER/ and /CWMS_WINTER_HOLDOUT/ | ELEV_TO_STOR_MOS.py experiment (script archived in dataprep). |
| //MOS/ELEV//1HOUR/CWMS_2011/ | Working record for the Jan 2011 CDB fill. Delete once you confirm the fill is folded into CWMS-CLEAN. |
| //MOS/HF//IR-MONTH/* (IRGZZAZD, IRXZZAZD, DFGZZ ZD, ...) | CWMS download version-flag artifacts; nothing reads them. |
| //MOS/STOR-COUNT//1DAY/CWMS/ (the copy in obsData) | QC output formerly written into the input store; MOS_STOR_RECORD_COUNT now writes it to MOS_Cleaned.dss. Regenerable. |
| //MAYFIELD OLD OUTFLOW/FLOW//1HOUR/USGS 14238000/ | Explicitly "OLD". |
| //CASTLE_ROCK/FLOW-UNREG//1DAY/SIMPLE/ | Output of the removed simplified/basic scripts. |
| /CASTLEROCK_NWS/FLOW-UNREG//1DAY + 1YEAR/BASIC_MOSSYROCK_STORAGE.../ | Output of removed #BasicUnregCAS.py. |

## DELETE — probably, verify first

| Record family | Caution |
|---|---|
| //MAYFIELD OUTFLOW/FLOW + STAGE (15MIN and 1HOUR, F=USGS 14238000) | Legacy naming; only the retired Cowlitz.py module referenced it. Verify the ResSim watershed time-series mapping doesn't point at these before deleting. |
| /14243000/, /14238000/, /14226500/ short-A-part variants and .../FLOW//1HOUR/USGS-PST/ | Older naming and a PST-shifted variant; no script reads them. If the PST record was made for a timezone check, it's regenerable. |
| /DUMMY/MOS/ELEV/, /DUMMY/DUMMY/ELEV//1HOUR/DUMMY ELEV/, /DUMMY ELEV/ZERO/... | Only ZERO_DUMMY/ZERO/FLOW is referenced by Create_Ensembles, but the ResSim watershed itself may map dummy ELEV records. Verify in ResSim before deleting. |

## Other DSS files

- `output/MOS_Cleaned.dss` — all records current (holdout raw/processed,
  routed pair, routed diff, unreg, and now STOR-COUNT). Re-running the
  holdout script refreshes them; nothing to prune, but a periodic
  squeeze keeps the file small.
- `ref_in/CastleRockPeaks.dss`, `ref_in/results.dss` — removed from the
  repo 24 Jul 2026: they are final outputs of the archived study, not
  inputs to anything here. They live on in Cowlitz_FF_DataPrep and the
  zip. Delete your local copies under CAS_Unreg_FF/data/ too if you
  carried them over.
- `src/Cowlitz_Unreg/Cowlitz/Cowlitz.py`, `CowlitzConfig.xlsx`,
  `ref_in/AllData.pickle` — legacy ResSim-support module and its data;
  nothing imports Cowlitz.py anymore (only utilsDSS/utilsTime/utilsIO/
  HydrologicRouting are imported). Removal candidates pending your OK.

## Input/output convention (adopted)

`data/` and `ref_in/` hold SOURCE records only: observed downloads plus
the hand-cleaned store. `output/` and `ref_out/` hold everything a
script writes — intermediates consumed by later scripts (MOS_Cleaned.dss)
and final products alike. The one deliberate exception: Reg_Unreg
inflow scripts write derived records (calc/cleaned/volcor inflows,
locals) into the shared obsData.dss because ResSim reads them from
there; those are flagged KEEP above.
