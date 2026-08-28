# COWLITZ_INDEX — paste or fetch this first in a new chat

Compact context for the Cowlitz flow frequency work in repo
`JoeShmoJo/Cowlitz_FF`. Read this INSTEAD of crawling the tree; fetch
individual files only as needed. Last updated 26 Aug 2026.

Full run order for both projects -- every script and every ResSim run, in
sequence, with a flow chart: `RUN_ORDER.md` at the repo root.

## The three project folders

| Folder | Purpose | Status |
|---|---|---|
| `CAS_Unreg_FF/` | Unreg flow frequency at Castle Rock (USGS 14243000) via MOS holdout — no ResSim | ACTIVE |
| `CAS_Reg_Unreg/` | ResSim inputs, the regulated curve, and the extension downstream to the Coweeman confluence | ACTIVE. Two memos in `docs/` |
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
wy_record_ssp.csv audit table). Source rules: pre-1968 Peak = USGS
peak record, 1/3/5-day = USGS daily record; post-1968 Peak = hourly
unreg else dS_2day regression, One_day = hourly 1-day avg else daily
unreg one-day max, 3/5-day = unreg daily averages. NOTHING reads from
the Cowlitz_FF_DataPrep archive; CastleRock_USGS_peaks.csv lives in
CAS_Unreg_FF/data. QC: `MOS_STOR_RECORD_COUNT.py`. CAS_Unreg_FF runs FIRST; it
owns obsData.dss and the download script.

## CAS_Reg_Unreg method (current)

Regulated curve: the unregulated curve is carried through a
reg-vs-unreg transform derived from ResSim, not fitted to the observed
regulated record (which is neither homogeneous nor log-Pearson).
Adjusted regulated peaks + synthetic ensembles build the transform;
`#Unreg_Reg_Curve.py` applies it -> `regulated_frequency_inferred.csv`.

Downstream extension (Section 8 of the combined memo). Four locations,
Castle Rock gage to the Coweeman confluence:

    Q(p) = Q_reg(p) + Q_unreg(p) x (DA_site - 2229)/2229 x 0.80

    Castle Rock gage        2229 sq mi        --
    below Arkansas Creek    2278 sq mi      49 sq mi   2.2%
    below Ostrander Creek   2335 sq mi     106 sq mi   4.8%
    below Coweeman River    2476 sq mi     247 sq mi  11.1%

The local term scales off the UNREGULATED curve because the tributaries
enter below the projects and respond to the storm, not the release.
Pairing unreg(p) with reg(p) is NOT a coincidence assumption -- same
river, same event, and reg was derived from unreg by routing. Nothing
assumes a tributary is simultaneously at its own AEP; that is what the
superseded coincident-frequency scripts did.

Incremental area, not named-tributary area: 247 sq mi against 197.5 for
the three named basins, the difference being ungaged local drainage.

Plain area ratio, no unit-runoff uplift, on three independent lines --
USGS 14245000 annual peaks WY1950-1996 converge to 1.04x area share at
WY1996 (212,245 cfs unreg, 92% of the 1,000-yr); Ecology 26C075 gives
1.11x for its top bin and that is a LOWER bound because its rating
ceiling censors the largest events; PRISM basin precipitation ratio is
1.04 over 71 years, so the two basins get the same depth and the
common-event excess is a response difference that vanishes as both
saturate.

Lag factor 0.80: measured Coweeman flow at the regulated crest over its
own event peak is 0.806 (n=52) and 0.781 (n=19) in the two well-sampled
bins, and 0.80 is also the 2009 study's value. The >60k bin gives 0.413
but on n=7, its mechanism is lead time rather than magnitude, and the
three events where the Coweeman exceeded its rating are missing from it.

1,000-yr below the Coweeman: 209,902 cfs (local 20,468, +10.8% over the
gage). Uncertainty is the gage's own band translated by the local term;
no extra uncertainty is added for the local.

## Critical facts

- `CAS_Unreg_FF/data/obsData.dss` is the canonical LIVE observed store,
  shared across projects. `//MOS/ELEV//1HOUR/CWMS-CLEAN/` in it is
  HAND-EDITED in DSSVue, irreproducible; some series manually copied
  in. Never regenerate/overwrite casually. CAS_Reg_Unreg scripts read
  and write this file via `../../CAS_Unreg_FF/data/obsData.dss`.
- No record duplication across projects — cross-project reads instead.
- `Cowlitz_FF_DataPrep/data/obsData.dss` (committed) = frozen memo-01
  snapshot; do not confuse with the live store.
- Daily MOS elevation for the regression + daily mass balance =
  //MOS/ELEV//1DAY/USGS/ (separate daily record, NOT derived from
  hourly). Hourly CWMS-CLEAN is used ONLY in the holdout workflow,
  gated day-by-day by STOR-COUNT + overrides -- much of the post-1974
  hourly record is unusable at hourly resolution. Legacy
  IRVZZAZD_CLEANED daily record retired; deletion list in
  CAS_Unreg_FF/docs/DSS_RECORD_CLEANUP.md.
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
- Ecology gage exports are FIXED WIDTH and leave the value column BLANK
  when nothing is reported, putting the quality code alone on the line.
  A regex that makes the value mandatory reads code 254 ("rating table
  exceeded, data will not be reported") as a discharge of 254 cfs --
  the gage's highest flows become one of its lowest values, silently.
  The parser lives ONCE in `Modules/ecology_io.py`; nothing re-declares
  it. Codes 151 and 254 carry no usable value and are NaN. Code 160 is
  in no legend in any file yet holds 1,012 readings up to the record
  maximum.
- Coweeman gage 26C075 rating tops out near 3,400 cfs. Of the nine
  events in the >60k bin, ZERO have a peak from an in-rating gaged
  reading. Any Coweeman tail statistic is a lower bound.
- `MEMO_CAS_Combined_FlowFrequency_DRAFT.docx` is HAND-EDITED and is the
  source of truth for prose. `make_memo_combined.py` regenerates it from
  source and will destroy those edits. DRAFT2 is built by
  `make_memo_draft2_section8.py`, which copies DRAFT and inserts only
  Section 8 + Appendix F.
- User environment: Windows, network-restricted (SSL inspection),
  runs scripts locally; container/Claude cannot run pydsstools —
  verify logic with stubs/synthetic data instead.

## Status: CAS_Unreg_FF RECORD COMPLETE (24 Jul 2026)

WY1927-2026 in output/CAS_Unreg_SSP.dss (import *-MOS-HOLDOUT paths in
SSP); only WY1969-1973 absent. Peak 95 (42 pre-reg USGS / 30 dS_2day
regression / 23 holdout), One_day 94, 3/5-day 93 (WY2026 3/5-day await
year close). Record made monotonic backwards from 5-day: 10 raises
across 9 WYs, logged in diagnostics/record_monotonic_adjustments.csv. Adopted:
(REG-UNREG) = -0.869*dS_2day - 3,836 cfs, R2 0.871, SE ~5,530, n 17.
Season override = full regulated era WY1974-2025 (winter-only daily
record; recorded season assumed to capture maxima); WY2026 awaits
year close. Memo final in CAS_Unreg_FF/docs.

## Open tasks

1. WY1969-1973 SSP treatment (perception threshold vs exclusion).
2. Review low-coverage WYs in the regression fit set.
3. WY2026 close-out re-run (3/5-day).
4. Prune retired DSS records per CAS_Unreg_FF/docs/DSS_RECORD_CLEANUP.md,
   then squeeze obsData.dss.
5. CAS_Reg_Unreg: implement cleaned-elevation-derived MOS inflows for
   ResSim (future).
6. Upload the hand-chopped `ensemble_synthetic.dss` to close the
   provenance gap in RUN_ORDER.md (results are post-chop, the committed
   ensemble is not).
7. Confluence map: verify the three approximate seed coordinates and the
   delineated basin area on a machine with NLDI reachable, then tune
   `label_frac` for the three creek names.
8. Section 2 of the combined memo cites 2,238 sq mi above Castle Rock;
   Section 8 and the StreamStats delineation use 2,229. Reconciled in
   text, not in the numbers.

## Key file map

- CAS_Unreg_FF: `README.md` (method + run order), `src/*`, output
  `MOS_Cleaned.dss` paths incl.
  `//CASTLE ROCK/FLOW-UNREG//1HOUR/CAS+ROUTED-DIFF/`,
  `wy_peak_records.csv`, `peakdiff_storage_regressions.csv`,
  `unreg_peak_estimates.csv`, `wy_record_ssp.csv`,
  `CAS_Unreg_SSP.dss` (final SSP input), memo in `docs/`.
- CAS_Reg_Unreg: `README.md` (script + diagnostics use cases);
  `output/regulated_frequency_inferred.csv` (the regulated curve),
  `output/below_confluence_frequency.csv` + `output/freq_table_*.csv`
  (the four downstream locations), `docs/Downstream_Confluence_Notes.md`
  (method history, evidence, and the Ecology parsing bug),
  `docs/CDID3_Coincident_Frequency_Notes.md`,
  `docs/MEMO_CAS_Reg_Unreg_DRAFT.docx`,
  `docs/MEMO_CAS_Combined_FlowFrequency_DRAFT.docx` (hand-edited) and
  `..._DRAFT2.docx` (adds Section 8); external ResSim watershed paths
  under `C:\Projects\Cowlitz_Flow_Frequency\ResSim\...` (not in repo).
- Shared: `Modules/ecology_io.py` (WA Ecology parser + quality codes).
- Cowlitz_FF_DataPrep: `RUN_ORDER.md` (archive note, diagnostics
  index, memo pipeline), `docs/MEMO_01_...docx`.
