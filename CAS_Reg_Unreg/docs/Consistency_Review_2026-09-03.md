# Consistency review, 3 September 2026

Reviewed:

    MEMO_CAS_Combined_FlowFrequency_2026_03_09.docx
    CowlitzFlowFreq_DQC_Memo_2026_03_09.docx

No document was edited. Every table in the memo was checked value by value
against the CSV that produces it, every embedded figure was checked against the
live script output, and the prose numbers were checked against both.

Both documents are clean of tracked changes and open comments.

---

## A. Fix before this goes out

### A1. The 5 percent and 95 percent headings mean opposite things in the same document

Appendix B puts the HIGHER flow under "5% limit" and the LOWER flow under
"95% limit". Table 5-3 and Appendix E do the reverse. All 64 rows of Appendix B
are the reverse of Table 5-3 and Appendix E. Not some, all.

At the 1 percent AEP:

    Appendix B, Peak   5% limit 210,386      95% limit 145,504
    Table 5-3          5% 83,114             95% 172,352
    Appendix E, T E-1  5% 83,114             95% 172,352

Both conventions are real. HEC-SSP labels a confidence limit by the exceedance
probability OF THE LIMIT, so its "5 percent limit" is the upper flow, and
Appendix B follows HEC-SSP faithfully. Table 5-3 and Appendix E label by the
percentile of the flow, so their 5 percent is the lower flow.

Section 5.7 settles which one the memo means. It says the limits are
"a 5th or 95th percentile". Under that sentence, 5 percent is the lower flow,
so Appendix B is the one that contradicts the document's own stated convention.

This is the exact confusion the DQC reviewer asked to have removed. Their
comment was that it "invites confusion to call it lower 95% and upper 95%,
recommend just saying 5% and 95% like SSP does". The regulated tables were
fixed and Appendix B was not, so the document now contains both readings of the
same two words. A reviewer who checks one appendix against the other will find
the tables appear to disagree by a factor of two.

Either flip the Appendix B column order, or head those two columns
"Upper confidence limit" and "Lower confidence limit" and say once that HEC-SSP
numbers them by exceedance probability.

Response: Switch the column headers in Appendix B. Remove the word limit, so jsut 5% and 95%.

### A2. Section 5.7 describes a method the equations and the code do not use

The prose says:

> Those are converted to a standard deviation, combined with the transform
> standard deviation by root-sum-of-squares, and the single z = 1.645 for a
> 90 percent two-sided interval is applied once to the combined value

That describes combining the two sigmas in log units and applying z once at the
end. The equations directly below it, and `#Unreg_Reg_Curve.py`, do something
else. They apply z to each source separately, convert each to a flow offset in
cfs, and combine the two cfs offsets by root sum of squares. That is
`FREQ_TERM_MODE = "literal"`, which is what produced every number in Table 5-3
and Appendix E.

The two are not interchangeable:

| AEP | Table 5-3 lower | Prose method | | Table 5-3 upper | Prose method | |
|---|---|---|---|---|---|---|
| 10% | 57,545 | 59,350 | +3.1% | 91,370 | 88,758 | -2.9% |
| 2% | 73,092 | 77,574 | +6.1% | 143,827 | 138,472 | -3.7% |
| 1% | 83,114 | 89,167 | +7.3% | 172,352 | 164,800 | -4.4% |
| 0.5% | 98,137 | 105,382 | +7.4% | 211,031 | 202,069 | -4.2% |
| 0.1% | 147,271 | 156,184 | +6.1% | 330,843 | 321,499 | -2.8% |

The equations are right and the sentence is wrong, so this is a prose fix, not a
recomputation. It matters because that sentence is the one a reviewer would use
to reproduce the band.

Response:  If the equation is correct I think I can leave it and remove the description and let the equations do the talking. 

### A3. The Lower5 equation adds where it should subtract

    Upper95 = RegBest + sqrt((Unreg95 - UnregBest)^2 + (Reg95 - RegBest)^2)
    Lower5  = RegBest + sqrt((UnregBest - Unreg5)^2 + (RegBest - Reg5)^2)

Both are built properly in the equation editor, with real radicals and
superscripts, so the reviewer's comment on that is resolved. But the second one
carries a plus. As written the lower limit is above the best estimate. The code
is `reg_lower = reg_curve - delta_lo`, so the sign should be a minus. The
numbers in the tables are correct, only the printed equation is wrong.

Response: Change the sign to a minus.

### A4. Table 5-2 and Figure 5-4 were never re-run after the screening was tightened

Table 5-2 reports n = 44 for water years 1974 through 2024. Section 5.2 says 51
years were assessed and 10 were screened out, leaving 41.

`critical_duration_adjusted_fits.csv` is dated 14 August and still carries four
water years that the current screening rejects, and is missing one it accepts:

    included but now screened out   WY1974, WY1985, WY2006, WY2013
    accepted but absent             WY1977

Two of those four, WY1974 and WY2013, were removed under `reg_over_unreg`. That
is the screen Section 5.2 describes as catching a case that is
"physically possible only if release exceeds inflow during a large flood". So
the critical duration table rests in part on two pairs the memo elsewhere calls
impossible.

Re-fitting on the current screened set does not change the conclusion. The peak
still ranks first and the numbers improve slightly:

| Duration | Table 5-2, n=44 | Re-fit, current screen, n=40 |
|---|---|---|
| Peak (1-hr) | log R² 0.847, slope 0.602 | log R² 0.879, slope 0.621 |
| 1-Day | 0.830, 0.605 | 0.860, 0.622 |
| 3-Day | 0.783, 0.633 | 0.802, 0.645 |
| 5-Day | 0.734, 0.663 | 0.743, 0.677 |

So the adopted duration is safe. Re-run `#Critical_Duration_Adjusted.py` and
refresh Table 5-2, Figure 5-4, and the n in the caption.

I re-ran te script. I'm not getting the same stats as you put in here. Check that, and update teh memo.
---

## B. Should fix

### B1. Five figures in the memo do not match the live script output

Run `#Figure_Check.py` in this folder for the current list. As of this review:

    Figure 4-1   document holds figures_combined/regression_2day.png
    Figure C-1   document holds figures_combined/unreg_curve_peak.png
    Figure C-2   document holds figures_combined/unreg_curve_1day.png
    Figure C-3   document holds figures_combined/unreg_curve_3day.png
    Figure C-4   document holds figures_combined/unreg_curve_5day.png

The four Appendix C figures predate the axis change. They still run the
probability axis from 0.999, and the current script starts it at 0.99.

Figure 4-1 is the opposite case and needs a decision rather than a re-paste.
The version in the memo is correct, since its R² of 0.871, standard error of
5,528 and slope of -0.869 all match Table 4-2 and
`peakdiff_storage_regressions.csv` exactly. But no script now produces it. The
current `PeakDiff_Storage_Regression.py` writes a visibly different figure at a
different size, and `figures_combined/regression_2day.png` is a hand renamed
copy of something else. The memo figure is right and unreproducible, which is
the worse half of this problem.

Response: This stuff doesn't matter I don't think. If the figures with the truncated x axis are in teh folder, add replace them. If not, just leave what's there. If it's not wrong, it's not worth chancing introducing an error.

### B2. The Figure 5-5 caption no longer describes the figure

Caption: "The twelve observed source storms, unscaled and aligned on each
storm's own peak." The figure now plots days into the member window, 0 to 30,
with no peak alignment. Drop the last clause.

Response: Drop the last clause.

### B3. The title inside Figure 5-5 contradicts its caption

The image is titled "Synthetic source storms, unscaled". These are the observed
hydrographs the synthetics are built from, which is what the caption below
correctly says. Change the title in `#Create_Synthetic_Ensembles.py`.

Response: Synthetic source storms are the source storms used to build synthetics. If this is confusing, rename to "Source storms, unscaled"

### B4. Table 7-2 does not match the current timing output

    memo      20-40k n=52 m=0.81   40-60k n=19 m=0.78   >60k n=7 m=0.41   all n=78 m=0.79
    current   20-40k n=51 m=0.81   40-60k n=18 m=0.76   >60k n=7 m=0.42   all n=76 m=0.79

Two events have left the record since the table was built. Section 7.3 also
says 78 events in the text. The adopted 0.80 factor is unaffected, since the
all-events median is 0.79 either way.

Response: fix.

### B5 to B8. The DQC memo

- The "Conclusions/Recommendations" heading has nothing under it.
- The section references are to the old numbering. The comments cite Section
  8.2, Section 8.4, Section 5.6, Appendix F and Figure 5-3. Those are now
  Section 7.2, Section 7.4, Section 5.7, Appendix E and Figure 5-7. Quoting the
  draft as reviewed is fair, but a reader checking a resolution will go looking
  for a Section 8 that does not exist. One added "now Section 7.2" per comment
  would fix it.
- It cites `MEMO_CAS_Unreg_FF_DQC_Comments.docx`, which is nowhere in the repo,
  and `MEMO_CAS_Combined_FlowFrequency_26Aug2026_DQC_comments.docx`, which has
  been moved to `docs/Archive/`.
- Two comment threads end on an Evaluation with no Backcheck, the Section 5.1
  detail request and the uncertainty band clipping. Both were acted on. For a
  final they need the reviewer's close out.

Response: Ignore
---

## C. Editorial

- The filenames read `2026_03_09` while both documents are dated 3 September
  2026, which is `09/03/2026` in the memo header. Anyone sorting by filename
  will read them as 9 March.
  
  Response: I will fix
  
- Four semicolons remain in body prose, in Section 4.8, the Figure 5-6 caption,
  the Appendix A note and the Appendix E note. One em dash is used as
  punctuation in Section 5, "flood risk management—it is also operated".
  
  Response: Leave it
  
- Section 5.1: "shown in table 5-1 ." Lowercase table, and a space before the
  period.
  
  # Response: Fix it
  
- Section 5.1: "but some it is likely the project would increase releases" is
  garbled.
  
  # Response: Remove "some"
  
- Section 5.1: "lead to reasonable validation" should be "led".
  # Response: Fix it
	
- Section 4.7: "Ten values across nine water years were raised ," has a space
  before the comma.
  
  # Response: Fix it
  
- Section 5.2 mixes number styles in three consecutive sentences: "Thirty-six",
  "Ten", then "5 years received no adjustment, 4 where".
  
  # Response: Leave it
  

---

## D. Checked and correct

Worth recording, because it is most of the document.

- **Table 5-3**, all 35 values against `regulated_frequency_inferred.csv`,
  including the reduction percentages and both band limits.
- **Table 7-3**, all 24 values against `below_confluence_frequency.csv`.
- **Appendix A**, all 95 water years and all four durations against
  `wy_record_ssp.csv`.
- **Appendix D**, all 51 water years against `adjusted_peaks.csv`.
- **Appendix E**, all four tables, 16 ordinates each, against the
  `freq_table_*.csv` set.
- **Appendix B**, the computed and expected columns, all 64 rows. Only the two
  limit columns are affected by A1.
- **Table 4-2** against `peakdiff_storage_regressions.csv`, and n = 17.
- **Table 4-3**, the source counts. 42 plus 30 plus 23 is 95 for the peak, and
  the durations reconcile.
- **Table 4-4** against `CAS_Unreg_distribution_parameters.csv`, every
  parameter and every count.
- **Section 4.1**, the 23 holdout years listed by water year count to 23.
- **Section 5.2**, the screening breakdown. Seven timing mismatches, two
  regulated above unregulated, and WY1980, matching
  `adjusted_peaks_screened_out.csv` exactly. 36 plus 10 plus 5 is 51.
- **Section 5.4**, the synthetic regulated peak range. The file gives 90,458 to
  266,744 against the memo's 90,000 to 267,000.
- **Section 5.6**, the reduction narrative. The maximum is 34.4 percent at the
  2 percent AEP and it falls to 16.3 percent at 0.1 percent.
- **Section 5.6**, the Figure 5-8 plotting position paragraph now matches
  `ADJUSTED_PP_BASIS = "ranked_unreg_positions"` in the script.
- **Section 7.1 and 7.2**, the drainage areas and incremental areas in
  Table 7-1, and the 0.80 timing factor throughout.
- **Figures 4-2, 5-4 through 5-8, 6-1, 6-2, 7-1 and 7-2** all match their live
  script output.

---

## E. One thing this review could not settle

The unregulated record changes estimation method at the same water year it
changes level. WY1927 to WY1968 are observed USGS pre regulation peaks with a
median of 53,850 cfs. WY1974 to WY2026 are reconstructed, either
`dS2day_regression` or `hourly_holdout`, with a median of 72,455 cfs.
Regulation began in 1968 and the record has a gap from WY1969 to WY1973, so
there are no overlap years and the two methods cannot be compared against each
other directly.

The step may be a real climate signal, since the Pacific Decadal Oscillation
shifted in the mid 1970s. It may be the reconstruction running high. The
frequency curve is fitted across the join either way, so if any part of the
step is method bias it propagates into every result in the memo. This is not a
finding, it is a question worth an answer before the next review.

# Response: Ignore it
