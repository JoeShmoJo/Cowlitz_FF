# Handoff: consistency review of the Cowlitz combined flow frequency memo

Written for a reviewing agent joining cold. Read this before touching anything.
Everything below was verified against the repo, not assumed.

Repo: `JoeShmoJo/Cowlitz_FF`. Working branch for this effort is
`claude/synthetic-scaling-multi-peak-x0q25a`.

Read `CLAUDE.md` at the repo root first. It carries the DSS, pydsstools, and
pandas traps and it is not optional. This file does not repeat it.

---

## 1. What the study produces

Peak flow frequency for the Cowlitz River at Castle Rock, unregulated and
regulated, plus an extension downstream to the Coweeman confluence. It supports
a levee system evaluation and it will be reviewed outside the district, so
internal consistency between the memo, the tables, the figures and the scripts
is the deliverable, not a nicety.

The chain, in the order the numbers are actually created:

```
observed and reconstructed record        CAS_Unreg_FF/output/wy_record_ssp.csv
  -> HEC-SSP Bulletin 17C analysis
  -> unregulated frequency curve         CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv
  -> LOESS unregulated to regulated transform
  -> regulated frequency curve           CAS_Reg_Unreg/output/regulated_frequency_inferred.csv
  -> drainage area and timing scaling
  -> downstream locations                CAS_Reg_Unreg/output/below_confluence_frequency.csv
```

Every regulated ordinate is `transform(unreg_expected(p))`. The regulated curve
is not fitted to regulated data and has no probability axis of its own. That one
fact resolves several apparent inconsistencies, so keep it in mind.

## 2. Documents under review

Both live in `CAS_Reg_Unreg/docs/`.

| File | What it is |
|---|---|
| `MEMO_CAS_Combined_FlowFrequency_comments_2026-09-03.docx` | The memo. Carries the author's replies to the DQC reviewer. Four reviewer comments remain open in the file. |
| `CowlitzFlowFreq_DQC_Memo_comments_2026-09-03.docx` | The DQC review memo with replies. |

**Do not edit either document.** The author is editing them by hand and pastes
figures in themselves. A review finding is a written finding, not an edit.

The memo structure as it stands:

```
1 Purpose                        5 Regulated Flow Frequency at Castle Rock
2 Basin Description                5.1 HEC-ResSim Flood Operation Calibration
  2.1 Previous Studies             5.2 Adjusted Regulated Peak Record
3 Data Used                        5.3 Critical Duration
4 Unregulated Flow Frequency       5.4 Synthetic Flood Ensembles
  4.1 Hourly Unregulated Record    5.5 Unregulated-Regulated Transform
  4.2 Peak by Storage Regression   5.6 Regulated Frequency Curve
  4.3 Durations, Daily Balance     5.7 Regulated Frequency Uncertainty
  4.4 Pre-Regulation Period      6 Comparison with the 2009 Restudy
  4.5 Record Extension           7 Extension to the Coweeman Confluence
  4.6 Assembled Record           Appendix A  Unregulated record by water year
  4.7 Peak Duration Consistency  Appendix B  Unregulated ordinates
  4.8 Flow Frequency Analysis    Appendix C  Unregulated curves by duration
                                 Appendix D  Adjusted regulated record by WY
                                 Appendix E  Regulated ordinates by location
```

## 3. Figure and table provenance

Figures are staged in `CAS_Reg_Unreg/docs/figures_combined/` and pasted into the
memo by hand. The staging folder is therefore **not** guaranteed to match what is
embedded in the docx, and right now it does not.

Verified by md5 against `word/media/` in the current docx:

| Staged figure | In the docx? |
|---|---|
| `below_confluence_frequency.png` | yes, image16 |
| `critical_duration.png` | yes, image8 |
| `reg_vs_2009.png` | yes, image14 |
| `regression_2day.png` | yes, image3 |
| `synthetic_events_overlay.png` | yes, image9 |
| `unreg_adopted_curves.png` | yes, image4 |
| `unreg_curve_peak/1day/3day/5day.png` | yes, image17 to image20 |
| `unreg_vs_2009.png` | yes, image13 |
| `basin.png` | no, the docx carries a jpeg |
| `coweeman_peak_ratio.png` | **no, docx is older** |
| `final_uncertainty.png` | **no** |
| `synthetic_event_single.png` | **no, docx is older** |
| `synthetic_events.png` | **no** |
| `transform_scatter.png` | **no** |
| `unreg_reg_final_uncertainty.png` | **no, docx is older** |

The four marked "docx is older" are figures regenerated in the last few commits
that the author has not yet pasted in. That is expected and is being tracked. It
is also the single largest consistency hazard in the package, so verify any
number you read off a memo figure against the current PNG before calling it a
finding.

Memo figure numbers map to files as follows. The numbering shifted by three when
Figures 5-1 to 5-3 were added, so older commit messages are off by that amount.

| Memo | File | Script |
|---|---|---|
| 4-1 | `regression_2day.png` | `CAS_Unreg_FF/src/PeakDiff_Storage_Regression.py` |
| 4-2 | `unreg_adopted_curves.png` | `CAS_Unreg_FF/src/Frequency_Curves_And_Table.py` |
| 5-1 to 5-3 | ResSim validation, author supplied | not scripted here |
| 5-4 | `critical_duration.png` | `#Critical_Duration_Adjusted.py` |
| 5-5 | `synthetic_events_overlay.png` | `#Create_Synthetic_Ensembles.py` |
| 5-6 | `synthetic_event_single.png` | `#Create_Synthetic_Ensembles.py` |
| 5-7 | `transform_convergence.png` | `#Transform_Convergence_Figure.py` |
| 5-8 | `unreg_reg_final_uncertainty.png` | `#Unreg_Reg_Curve.py` |
| 6-1 | `unreg_vs_2009.png` | `Frequency_Curves_And_Table.py` |
| 6-2 | `reg_vs_2009.png` | `#Unreg_Reg_Curve.py` |
| 7-1 | `coweeman_peak_ratio.png` | `#Coweeman_HistoricPeakRatio.py` |
| 7-2 | `below_confluence_frequency.png` | `#BelowConfluence_FlowFrequency.py` |

Tables to check against their source CSVs:

| Memo table | Source |
|---|---|
| 4-3 assembled record | `CAS_Unreg_FF/output/wy_record_ssp.csv` |
| 4-4 distribution parameters | `CAS_Unreg_FF/output/CAS_Unreg_distribution_parameters.csv` |
| 5-1 adjusted vs unregulated by duration | `CAS_Reg_Unreg/output/critical_duration_adjusted_fits.csv` |
| 5-2 adopted regulated curve | `CAS_Reg_Unreg/output/regulated_frequency_inferred.csv` |
| 7-3 regulated peak by location | `CAS_Reg_Unreg/output/below_confluence_frequency.csv` |
| Appendix B | `CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv` |
| Appendix D | `CAS_Reg_Unreg/output/adjusted_peaks.csv` |
| Appendix E | `CAS_Reg_Unreg/output/freq_table_*.csv` |

## 4. What the consistency review should actually check

In rough order of how likely it is to find something real.

1. **Numbers quoted in prose against the CSVs.** The memo states many values
   inline. Each one should be traceable to a file in `output/`. Check the flow
   values, the percentages, the counts of water years and the record spans.
2. **Figure axes and captions against what the script draws.** Several figures
   were re-cut recently and a caption can survive a change to the plot.
3. **Table 5-2 and Appendix E against `regulated_frequency_inferred.csv`.** The
   regulated curve was refit after the synthetic ResSim run was corrected. Confirm
   both carry the current ordinates.
4. **Internal arithmetic.** Reduction percentages, drainage area ratios, the
   incremental area logic in Section 7, and the 0.80 timing factor.
5. **Record counts.** 95 unregulated water years, 51 assessed regulated water
   years spanning WY1974 to WY2024, 41 usable after screening. Any other count in
   the text is a finding.
6. **Terminology.** Expected probability curve versus computed curve. The memo
   adopts the **expected** curve throughout. The computed curve may appear in
   tables for reference only. Mixing them was a real defect once already.
7. **Uncertainty language.** The adopted band is a **90 percent two sided**
   interval, z equals 1.645, so the limits are the 5th and 95th percentiles.
   Column names in the CSV are `reg_lower_90pct_cfs` and `reg_upper_90pct_cfs`.
   Any surviving reference to a 95 percent band, or to z equals 1.960, is a
   finding. The combination follows EM 1110-2-1619 Section 4-4.c(1): each source
   is reduced to its own sigma in log10 units, combined by root sum of squares,
   and a single z is applied once at the end.
8. **The author's stylistic constraint.** No hyphens and no semicolons in memo
   prose. More sentences instead. Flag violations, do not fix them.

## 5. Things already settled. Do not relitigate.

- **Figure 5-8 plotting positions.** Settled on 3 Sep 2026 after two
  false starts. The DQC reviewer asked for the adjusted regulated peaks to be
  ranked so they increase continuously. Adopted:
  `ADJUSTED_PP_BASIS = "ranked_unreg_positions"` in `#Unreg_Reg_Curve.py`. The
  41 usable peaks are sorted largest first and given, in the same order, the
  plotting positions their own water years hold in the 95 year unregulated
  record. Monotone, on the curve's own probability axis, and centred on the
  curve (median observed over curve 0.99, 17 of 41 above). No value changes.
  The unregulated record stays as HEC-SSP draws it, all 95 years ranked once
  (`FINAL_UNREG_POINTS_BASIS = "full"`). Two things NOT to do again: ranking
  the regulated peaks inside the 51 year window puts the cloud above the curve
  (1.05, 33 of 41 above) because WY1974 to WY2024 holds nine of the ten largest
  unregulated years; and ranking the window's unregulated peaks alongside as a
  control makes the unregulated curve look like a poor fit to its own data,
  which invites a question the memo does not need. Both are kept only as
  diagnostics: `src/#Fig58_Placement_Options.py` draws all three placements.
  The memo carries one explanation paragraph before the figure that the author
  may delete, so nothing else references it.
- **Figure 5-7 lower band edge.** Beyond the edge of the fitted data
  (279,538 cfs) the lower edge of the transform scatter band closes its gap to
  the 1:1 line at the same average rate the drawn limb closes its own, and
  lands on 1:1 at about 393,000 cfs, off the plot. `LOWER_BAND_CLOSES_ON_1TO1`
  in `#Transform_Convergence_Figure.py`. The earlier parallel hold past the
  325,000 cfs convergence left a dogleg that read as the band widening.
  Figure only, no table changes.
- **The synthetic ensemble is hand edited after it is built.** See `CLAUDE.md`.
  `#Create_Synthetic_Ensembles.py` has `PLOTS_ONLY = True` deliberately. Do not
  set it to False. Running it in write mode destroys the author's hand chop of the
  Dec1977 and Nov1986 members and nothing downstream notices.
- **The transform figure solid to dashed split.** Cosmetic only. The
  `extrapolated` flag was always correct.
- **A scripted alternative to the hand chop** was tried in Aug 2026 and reverted.
  Do not re-propose it. `CLAUDE.md` has the reasoning.

## 6. Known open items, already on the author's list

These are flagged so you do not spend effort rediscovering them.

1. The 1906 historical event rationale, DQC memo comment 2. Needs the author.
2. Comment 252 anchors to an empty paragraph. Needs the author.
3. Comment 228 equations are correct as plain text with a variable list but still
   need rebuilding in the Word equation editor.
4. Comment 239 asks for a passage to move to Section 5.3. Deliberately not done,
   flagged as an editorial judgement.
5. `FREQ_TERM_MODE`. Both the literal and the `transform_curve` forms are computed
   on every run, the alternative sitting in the `_alt` columns of
   `regulated_frequency_inferred.csv`. Band width is unchanged at the median and
   ranges from 82 to 123 percent of current across the ordinates. Not adopted.
6. Figure 5-7's drawn convergence disagrees with Appendix E's tabulated
   extrapolated ordinates by up to 5.0 percent. At the 2,000 year the table says
   256,186 and the figure draws 268,935. Not reconciled.
7. `#MOS_Special_Release_MinFloodPool.py` still points at
   `CAS_Reg_Unreg/data/ObsData_RegUnreg.dss`, which was deleted. This is
   deliberate and is explained in `CLAUDE.md`.

## 7. A genuine open question found while preparing this handoff

Not yet raised with the author, and worth a look.

The unregulated record changes estimation method at exactly the same place it
changes level. WY1927 to WY1968 are observed USGS pre regulation peaks with a
median of 53,850 cfs. WY1974 to WY2026 are reconstructed, either
`dS2day_regression` or `hourly_holdout`, with a median of 72,455 cfs. There are
no overlap years, because regulation began in 1968 and the record has a gap from
WY1969 to WY1973, so the two methods cannot be compared directly against each
other.

That step could be a real climate signal, since the Pacific Decadal Oscillation
shifted in the mid 1970s. It could equally be a reconstruction running high. It
could be both. The frequency curve is fitted across the join, so if any part of
the step is method bias it propagates into every result in the memo.

Reproduce with:

```
python3 -c "
import pandas as pd, numpy as np
u = pd.read_csv('CAS_Unreg_FF/output/wy_record_ssp.csv')
print(u.groupby('Peak_Source').Peak.agg(['count','median','max']))
"
```

## 8. Ground rules

- Do not edit either .docx. Report findings in text.
- Do not run `#Create_Synthetic_Ensembles.py` with `PLOTS_ONLY = False`.
- Do not re-run `CAS_Unreg_FF/src/Frequency_Curves_And_Table.py` casually. It has
  a truncation guard because it once cut the frequency table from 16 ordinates to
  12. The guard refuses to shrink the table and prints which duration was short.
- `CAS_Unreg_FF/data/obsData.dss` is the shared store for both projects and
  several scripts write into it. Back it up before running those.
- Scripts `os.chdir` to their own folder and use relative paths from there, so
  run them from `src/`.
- Shared modules live in `/Modules` at the repo root. Nothing is duplicated.
- Anything that touches pydsstools must work on the author's Windows build too.
  Prefer `getattr` with fallbacks over a direct API call. `CLAUDE.md` lists the
  three API differences that have each already cost a round trip.
