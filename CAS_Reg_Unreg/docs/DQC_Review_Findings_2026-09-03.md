# Review of the 3 Sep 2026 DQC round

Findings only. Neither .docx was edited. Every number below was checked against
the CSVs in `output/`, not read off a figure.

Documents reviewed:

- `MEMO_CAS_Combined_FlowFrequency_comments_2026-09-03.docx` (the memo, four
  reviewer comments still open: C5, C18, C21, C31)
- `CowlitzFlowFreq_DQC_Memo_comments_2026-09-03.docx` (the DQC memo with
  evaluations and backchecks)

---

## 1. Figure 5-8 and comment C18: is the reviewer misguided?

The comment: "The regulated points should not be plotted this way. They should
be re-ordered so they continuously increase instead of plotting the regulated
peak at the same AEP as the unregulated peak for a given year."

**Short answer: no, the reviewer is not misguided, and neither is the author's
worry about the ranked version. Both are right about different things, and
there is a presentation that satisfies both.**

The reviewer is asking for the standard convention. A frequency curve is
checked against its record by ranking the record and plotting it at plotting
positions. That is what Appendix C does for the unregulated curves and it is
what every reader of this memo will expect for the regulated curve. Plotting a
regulated peak at its unregulated year's AEP is a legitimate diagnostic of the
transform scatter, but it is not a frequency plot, and it needs a paragraph of
explanation that the ranked plot does not.

The author's worry is also real. When the adjusted peaks are ranked (commit
`129acc5`) the cloud sits above the regulated curve and looks as though it could
not have produced it. Measured at median plotting positions against
`regulated_frequency_inferred.csv`:

| Placement of the 41 adjusted peaks | median obs / curve | above the curve |
|---|---|---|
| ranked among 41 usable values | 1.13 | 40 of 41 |
| ranked over the 51 assessed years | 1.05 | 33 of 41 |
| at each year's own unregulated AEP (current figure) | 1.03 | 22 of 41 |

**But the offset is a property of the sample window, not of the transform or
the curve.** The regulated record is WY1974 to WY2024. That window holds 9 of
the 10 largest unregulated years in the 95 year record, and its median
unregulated peak is 72,455 cfs against 61,703 cfs for the full record. Rank
anything from that window among itself and it will plot above any curve fitted
to the full record.

The proof is to do the same thing to the unregulated record, which nobody
disputes. Ranking the 51 unregulated peaks of WY1974 to WY2024 among themselves
and plotting them against the adopted unregulated curve:

| Record ranked inside WY1974 to WY2024 | median obs / curve | above the curve |
|---|---|---|
| unregulated, n = 51 | 1.16 | 43 of 51 |
| adjusted regulated, 41 of 51 | 1.05 | 33 of 41 |
| unregulated, full record for contrast, n = 95 | 1.01 | 49 of 95 |

The unregulated subset sits further above its own curve than the regulated
points sit above theirs. Nobody would conclude from that that the unregulated
curve is wrong. The regulated cloud is in fact the better behaved of the two.

`src/#Fig58_Window_Ranked_Check.py` draws this. Output is
`output/diagnostics/fig58_window_ranked_check.png`. It reads only published
CSVs and does not touch Figure 5-8.

**Recommendation.** Give the reviewer the ranked plot, because that is the
convention, and make the window effect visible on the same figure so that the
ranked cloud stops looking like a contradiction:

1. Rank the adjusted peaks over the 51 assessed years (n = 51, the ten screened
   years unranked), as `129acc5` did.
2. On the same figure, rank the WY1974 to WY2024 unregulated peaks among
   themselves and plot them in the unregulated colour. Keep the full 95 year
   unregulated record as a faint third series, or drop it.
3. One sentence in the caption or the text: both records are ranked within the
   51 year regulated period, and both plot above their curves because that
   period holds nine of the ten largest unregulated years in the 95 year
   record. The offset is the sample window, not the transform.

This does what the reviewer asked, answers the "cannot produce the curve"
objection with a control the reviewer already accepts, and drops the paragraph
that the own year placement would otherwise need. The own year placement can
stay as a diagnostic in `output/diagnostics/` if wanted, but it should not be
the memo figure.

Two things worth saying in the reply to the reviewer, briefly:

- The regulated curve has no probability axis of its own. Each ordinate is the
  transform of the unregulated ordinate at that AEP, so the ranked regulated
  points are a check, not the basis of the curve.
- The window effect is confounded with a method change. WY1927 to WY1968 are
  observed pre regulation peaks, WY1974 on are reconstructed. See item 7 of
  `HANDOFF_Consistency_Review.md`. That does not change the plotting position
  argument, since both the curve and the control use the same record.

Handoff section 5 records the own year basis as settled. This revisits it at
the author's request. `ADJUSTED_PP_BASIS = "assessed_span"` in
`#Unreg_Reg_Curve.py` already produces step 1. Step 2 needs the WY1974 to
WY2024 unregulated subset added as a series, which the diagnostic script shows
how to do.

---

## 2. The other three open memo comments

**C5, Figure 5-6.** The staged `synthetic_event_single.png` already has the
title "Scaling one source storm: Feb1996" and an x axis labelled "Days", and
the strain term is gone from the legend. Two things the reviewer may still ask
about: the title also carries "212,246 cfs unreg, WCM_RC ratio 0.62 -- largest
on record", and the legend entry for the fourth member is "beyond". Suggest
trimming the title to the storm name and labelling the fourth member "20%
above 500yr" to match Section 5.4. Not yet pasted into the docx.

**C21, Section 5.7.** The reviewer is right. The sentence says Appendix E lists
"the two terms that produced its uncertainty band" and Table E-1 has only
Unregulated, Regulated, 5% and 95%. `regulated_frequency_inferred.csv` carries
both terms (`sigma_freq_lo_dex`, `sigma_freq_hi_dex`, `sigma_transform_lo_dex`,
`sigma_transform_hi_dex`, and the cfs half widths `freq_term_*_cfs` and
`transform_term_*_cfs`). Either add two columns to Table E-1 or delete the
clause. Adding the columns is the better answer to the DQC comment that
started this, which asked how the bounds were drawn.

**C31, Section 7.4.** "Appendix F" should be "Appendix E". Same paragraph says
"table 8-3", which should be "Table 7-3".

---

## 3. DQC memo (CowlitzFlowFreq_DQC_Memo)

1. **The Review Description paragraph is boilerplate from another review.** It
   says the memo "provides the DQC comments on the Chapter 4 (Hydrologic
   Hazards) document for the Lost Creek Lake dam safety risk assessment". This
   will be the first thing an outside reviewer reads.
2. The header date is 03Aug2026 but the memo contains comments drafted 26 Aug
   and backchecks from 2 Sep. Update the date or add a revision line.
3. **Conclusions/Recommendations is empty.**
4. Comments cite the section numbers of the draft they were written against
   (Section 8.2, 8.4, Figure 5-3, Section 5.6, Appendix F). That is fine as a
   record, but each Evaluation could name the current location (7.2, 7.4,
   Figure 5-7, Section 5.7, Appendix E) so a later reader can find it.
5. Backcheck on the 5.1 comment asks for the numbers used for minimum release,
   maximum evacuation release, and the downstream control flow target, a
   statement that no special curves are published, and the assumed time of
   recession. None of that is in Section 5.1 yet. Open.
6. Backcheck on the uncertainty comment asks for the upper bound to be clipped
   at 1:1. The staged `transform_scatter.png` does this. The docx still carries
   the older figure where the band crosses 1:1 near 250,000 cfs, so the
   backcheck will look unanswered until the figure is pasted.
7. Backcheck on the synthetic events, the 1906 event, the Coweeman scatter
   plot, the Section 8.4 equation, and the 90 percent interval all concur.
   The memo text is consistent with each of those resolutions.

---

## 4. Memo consistency findings, in order of how much they matter

Checked against `wy_record_ssp.csv`, `adjusted_peaks.csv`,
`critical_duration_adjusted_fits.csv`, `regulated_frequency_inferred.csv`,
`below_confluence_frequency.csv`, `synthetic_results.csv`,
`CAS_Unreg_distribution_parameters.csv`, `CAS_Unreg_frequency_table.csv`,
`peakdiff_storage_regressions.csv`, and `coweeman_lagfactor_events.csv`.

1. **Section 5.2 counts do not add up.** "Thirty-eight received an adjustment
   ... Ten were screened out ... 5 years received no adjustment" sums to 53
   for 51 assessed years. From `adjusted_peaks.csv`: 41 usable, of which 36
   were adjusted and 5 were not (1978, 1983, 1998, 2020 with a zero
   difference, and 1976 with a negative one), plus 10 screened. The 38 most
   likely counted screened years that also had a positive difference. The
   median of the 36 adjustments is 10,300 cfs, not 11,000. The range 263 to
   24,013 cfs is stated as "250-24,000", which is close enough but 260 would
   be exact.
2. **Table 5-1 uses 44 pairs while the transform uses 41.** The critical
   duration dataset keeps WY1974, 1985, 2006 and 2013, all of which the
   adjusted record screens out. Either rerun `#Critical_Duration_Adjusted.py`
   on the screened record or add a sentence saying the duration comparison
   was made before screening. A reviewer who cross reads Table 5-1 and
   Figure 5-7 will ask.
3. **Table 4-4 is stale for the 3 day and 5 day columns and for the peak
   equivalent record length.** Against `CAS_Unreg_distribution_parameters.csv`:

   | Row | Memo 3-Day | CSV | Memo 5-Day | CSV |
   |---|---|---|---|---|
   | Mean of logs | 4.660 | 4.663 | 4.593 | 4.596 |
   | Std dev of logs | 0.192 | 0.194 | 0.179 | 0.181 |
   | Station skew | -0.076 | -0.079 | -0.051 | -0.054 |
   | Weighted skew | -0.007 | -0.009 | -0.029 | -0.031 |

   Peak equivalent record length is 103.7 in the memo and 100.7 in the CSV.
   The prose in 4.8 ("-0.158 to -0.054", "100-112 years") matches the CSV, and
   Appendix B matches `CAS_Unreg_frequency_table.csv`, so the table is the
   odd one out. Probably from the run before WY2026 was added.
4. **Section 5.4: "regulated peaks of 81,000-267,000 cfs".**
   `synthetic_results.csv` spans 90,458 to 266,744. The 81,000 is from an
   earlier run.
5. **Cross references that went stale when 5.1 was inserted.**
   - Section 5 intro: "(5.1)" adjusted peaks is now 5.2, "(5.2)" pairing is
     5.3, "(5.3)" synthetic floods is 5.4, "(5.4-5.5)" is 5.5 to 5.6.
   - Section 5.5: "Section 5.8 combines the regulation uncertainty" should be
     5.7.
   - Table 3-1 purpose column cites 5.2, 5.3 and 5.4 for the hourly and daily
     unregulated flow work that lives in 4.1 to 4.3.
   - Section 5.1: "Figures 5-1, 5-1, and 5-3".
   - Section 7.4: "table 8-3" and "Appendix F" (C31).
6. **Section 6.2 cites something that is not there.** "as noted in section
   5.5, has a shape that is likely to cause the reservoir to fill and spill"
   refers to the 1933 event. Neither 5.4 nor 5.5 mentions 1933. Either add the
   note where the source storms are described or drop the reference.
7. **Section 7.3: "average scaling factor of 0.79".** From
   `coweeman_lagfactor_events.csv` the median over 78 events is 0.79 and the
   mean is 0.75. Table 7-2 correctly says median. Change "average" to
   "median".
8. **Drainage area at the gage.** Section 7, Table 7-1 and the script use
   2,238 square miles. Table E-1's caption says 2,229. Section 2 gives 2,480
   for the total basin and Table 7-1 gives 2,476 below the Coweeman, which are
   different points and are fine if the text says so.
9. Section 4.1: "2008-2026" as holdout years includes 2013, which is a
   regression year. The count of 23 is right. Say "2008 through 2026 except
   2013".
10. Section 4.7: "eight of them regression-estimated peaks lifted to their
    1-day value". Eight peaks were lifted, but WY2019 is a holdout peak, so
    seven regression and one holdout.
11. Section 5.7 equation lines are placeholders (Upper95 =, Lower5 =, with a
    variable list below). Handoff item 3, still open.
12. Style, since the handoff flags it: an em dash in the Section 5 opening
    paragraph ("flood risk management—it is also") and in the Appendix E
    introduction ("local contribution — see Section 7.4"). Semicolons in the
    4.8 prose, the Appendix A note, the Figure 5-6 caption, and several Table
    3-1 and Table 4-3 cells.

Figures still to be pasted (staged PNG differs from the docx): 5-6, 5-7, 5-8,
and 7-1. The docx 5-8 is the own year version, so whichever of the section 1
options is adopted it needs re pasting.

---

## 5. Verified and consistent

- Table 5-2 and Table E-1 against `regulated_frequency_inferred.csv`, every
  ordinate.
- Table 7-3 and Tables E-2 to E-4 against `below_confluence_frequency.csv`.
- Appendix B against `CAS_Unreg_frequency_table.csv`, all four durations.
- Appendix D against `adjusted_peaks.csv`.
- Table 4-2 against `peakdiff_storage_regressions.csv`, n = 17 pairs.
- Table 4-3 counts: 42 + 30 + 23 = 95 peaks, 41 + 30 + 23 = 94 one day, 41 +
  53 = 94 three and five day.
- Section 5.6 reductions: 34 percent at 2 percent AEP, 16 percent at 0.1
  percent AEP.
- Section 4.7 largest monotonic changes: 64, 19 and 19 percent in WY2001,
  1985 and 2006.
- Uncertainty language: 90 percent two sided, z = 1.645, 5th and 95th
  percentiles, throughout.
- Record counts: 95 unregulated years, 51 assessed, 41 usable.

---

## 6. Known and deliberately not reopened

From the handoff: the synthetic ensemble hand chop, the scripted alternative
to it, `FREQ_TERM_MODE`, the Figure 5-7 versus Appendix E extrapolated ordinate
gap of up to 5 percent (still open, 256,186 in the table against 268,935 drawn
at the 2,000 year), and the pre 1968 versus post 1974 level step in the
unregulated record. The last one is worth a sentence in Section 4.6 before an
outside reviewer finds it, whatever its cause.
