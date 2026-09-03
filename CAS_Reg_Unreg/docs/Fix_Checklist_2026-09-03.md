# Fix checklist for the 3 September documents

Ordered so nothing has to be done twice. Findings are in
`Consistency_Review_2026-09-03.md`.

Marked **[you]** where it has to be done on your machine, either because it is
a Word edit or because the script needs a file that is not in the repo.

---

## Step 1. Fix the uncertainty method sentence in Section 5.7   [you]

**Finding A2. Ten minutes. No numbers change.**

The sentence describes combining standard deviations and applying z once. The
equations below it, and the code, combine flow offsets in cfs. The equations
are right, so only the sentence is wrong.

Find this:

> HEC-SSP reports its own limits at 5 and 95 percent. Those are converted to a
> standard deviation, combined with the transform standard deviation by
> root-sum-of-squares, and the single z = 1.645 for a 90 percent two-sided
> interval is applied once to the combined value, so every limit reported here
> is a 5th or 95th percentile.

Replace with:

> HEC-SSP reports its own limits at 5 and 95 percent. Each source is expressed
> as a flow offset from its own best estimate at the 90 percent two sided
> level, using z = 1.645 in log space. The two offsets are then combined by
> root sum of squares and applied to the regulated best estimate, so every
> limit reported here is a 5th or 95th percentile.

That now matches both the equations and `FREQ_TERM_MODE = "literal"` in
`#Unreg_Reg_Curve.py`.

---

## Step 2. Fix the sign in the Lower5 equation   [you]

**Finding A3. Two minutes. No numbers change.**

In the equation editor, second equation in Section 5.7. Change the leading plus
to a minus:

    Lower5 = RegBest - sqrt((UnregBest - Unreg5)^2 + (RegBest - Reg5)^2)

Leave Upper95 as it is. The code is `reg_lower = reg_curve - delta_lo`, so the
minus is correct and every number already in the tables is right.

---

## Step 3. Decide the 5 and 95 percent convention, then fix Appendix B   [you]

**Finding A1. Thirty minutes. This is the important one.**

Appendix B currently puts the higher flow under "5% limit". Table 5-3 and
Appendix E put the lower flow under "5%". All 64 Appendix B rows are the
reverse of the other two.

**Recommended fix.** Make Appendix B match Table 5-3 and Appendix E, because
Section 5.7 defines the convention as percentiles of the flow and those are the
adopted results. Swap the two column contents in Appendix B so the LOWER flow
sits under "5% limit", then add one sentence to the Appendix B lead in:

> HEC-SSP numbers its confidence limits by the exceedance probability of the
> limit itself, so what it prints as the 0.05 limit is the upper flow. The
> columns here are ordered by percentile of flow instead, matching Table 5-3
> and Appendix E.

That keeps the reviewer's "just say 5 and 95 percent" and warns anyone who
compares against raw SSP output.

**Check when done.** At the 1 percent AEP, Peak duration, Appendix B should
read 145,504 under 5% and 210,386 under 95%. Right now it reads them the other
way round.

**The alternative**, if you would rather stay literally faithful to SSP, is to
flip Table 5-3 and Appendix E instead. Do not do that lightly. It contradicts
the "5th or 95th percentile" wording in Section 5.7 and it means re-labelling
three tables instead of one.

---

## Step 4. Re-run the critical duration analysis   [you, Windows only]

**Finding A4. Thirty minutes including the re-paste.**

`critical_duration_adjusted_fits.csv` is from 14 August and predates the
tightened screening. It carries four water years that are now screened out, two
of them under `reg_over_unreg`, and misses WY1977 which is now accepted.

    cd CAS_Reg_Unreg\src
    python "#Critical_Duration_Adjusted.py"

It has to run on your machine because it reads the ResSim simulation.dss under
`C:\Projects\2026_Cowlitz_Flow_Frequency\`, which is not in the repo.

Then update **Table 5-2** and re-paste **Figure 5-4**. Expected new values:

| Duration | now in Table 5-2 | after the re-run |
|---|---|---|
| Peak (1-hr) | n 44, R² 0.870, log R² 0.847, slope 0.602 | n 41, log R² near 0.879, slope near 0.621 |
| 1-Day | 0.853, 0.830, 0.605 | near 0.860, 0.622 |
| 3-Day | 0.792, 0.783, 0.633 | near 0.802, 0.645 |
| 5-Day | 0.733, 0.734, 0.663 | near 0.743, 0.677 |

Change the n in the table and anywhere the caption or text says 44. The adopted
duration does not change. The peak still ranks first and every R² improves.

**Then re-run the curve, which is quick and safe:**

    python "#Unreg_Reg_Curve.py"
    python "#BelowConfluence_FlowFrequency.py"

I tested the effect of adding WY1977 back. **Every regulated flow in Table 5-3
is unchanged at every AEP.** The only visible movement is the 2-year band row,
which goes from 39,578 to about 39,914 on the lower limit and 60,614 to about
60,913 on the upper. Every other band value moves by tens of cfs. Table 7-3 and
Appendix E do not move at all. So this does not cascade, and you do not need to
revisit Section 5.6 or Section 7.

Re-paste **Figure 5-8** and **Figure 5-7** after the re-run.

---

## Step 5. Fix the Figure 5-5 title in the script   [me or you, one line]

**Finding B3. Two minutes.**

In `CAS_Reg_Unreg/src/#Create_Synthetic_Ensembles.py`, `plot_events_overlay()`:

    ax.set_title("Synthetic source storms, unscaled", fontsize=11)

These are the observed hydrographs the synthetics are built from, which is what
your own caption says. Change it to:

    ax.set_title("Observed source storms, unscaled", fontsize=11)

Then re-run the script and re-paste Figure 5-5. `PLOTS_ONLY` is True, so it
will not touch `ensemble_synthetic.dss`. Leave it True.

I did not make this edit myself because you are hand editing that file and I do
not want to hand you another merge conflict.

---

## Step 6. Fix the Figure 5-5 caption   [you]

**Finding B2. One minute.**

Current:

> Figure 5-5. The twelve observed source storms, unscaled and aligned on each
> storm's own peak.

The figure now plots days into the member window with no peak alignment. Cut
the last clause:

> Figure 5-5. The twelve observed source storms, unscaled.

---

## Step 7. Re-paste the four Appendix C figures   [you]

**Finding B1. Ten minutes.**

They still run the probability axis from 0.999. The current script starts it at
0.99. No re-run needed, the current output is already correct on disk. Stage
and paste:

    cd CAS_Reg_Unreg\docs
    python "#Figure_Check.py"

Then paste these four from `figures_combined/`:

    Fig_C-1__CAS_Unreg_freq_Peak.png
    Fig_C-2__CAS_Unreg_freq_1-Day.png
    Fig_C-3__CAS_Unreg_freq_3-Day.png
    Fig_C-4__CAS_Unreg_freq_5-Day.png

Run `#Figure_Check.py` again afterwards. All four should report CURRENT.

---

## Step 8. Decide what to do about Figure 4-1   [you]

**Finding B1, the awkward half. Twenty minutes.**

The figure in the memo is **correct**. Its R² of 0.871, standard error of 5,528
and slope of -0.869 all match Table 4-2 and
`peakdiff_storage_regressions.csv` exactly. But no script now produces it.
`figures_combined/regression_2day.png` is a hand renamed copy of something, and
`PeakDiff_Storage_Regression.py` currently writes a different looking figure at
a different size.

    cd CAS_Unreg_FF\src
    python "PeakDiff_Storage_Regression.py"

Then open `CAS_Unreg_FF/output/diagnostics/peakdiff_storage_regression.png` and
compare it against Figure 4-1 in the memo. If the new one is acceptable, paste
it and you are done. If the memo version is the one you want, the script needs
editing to reproduce it, because right now a reviewer asking "where did this
figure come from" has no answer.

Nothing in the numbers is wrong either way. This is about reproducibility.

---

## Step 9. Refresh Table 7-2 and the event count   [you]

**Finding B4. Fifteen minutes.**

Current output of `#Coweeman_RegPeak_Timing.py`:

| Castle Rock unregulated peak | Events | Median ratio |
|---|---|---|
| 20,000 to 40,000 cfs | 51 | 0.81 |
| 40,000 to 60,000 cfs | 18 | 0.76 |
| Above 60,000 cfs | 7 | 0.42 |
| All events | 76 | 0.79 |

Section 7.3 also says 78 events in the text. Change it to 76. The adopted 0.80
factor does not move, since the all events median is 0.79 either way.

---

## Step 10. Close out the DQC memo   [you]

**Findings B5 to B8. Twenty minutes.**

1. Write the Conclusions and Recommendations section. It is currently an empty
   heading.
2. Add the current section number to each comment that cites the old
   numbering. Section 8.2 is now 7.2, Section 8.4 is now 7.4, Section 5.6 is
   now 5.7, Appendix F is now Appendix E, and Figure 5-3 is now Figure 5-7.
3. Fix the two file references. `MEMO_CAS_Unreg_FF_DQC_Comments.docx` is not
   anywhere in the repo, and
   `MEMO_CAS_Combined_FlowFrequency_26Aug2026_DQC_comments.docx` has moved to
   `docs/Archive/`.
4. Get Ryan's backcheck on the two threads that end on an Evaluation with no
   Backcheck, the Section 5.1 detail request and the uncertainty band clipping.
   Both were acted on.

---

## Step 11. Editorial sweep   [you]

**Section C. Fifteen minutes. Find and fix.**

| Where | Now | Should be |
|---|---|---|
| Section 5.1 | "shown in table 5-1 ." | "shown in Table 5-1." |
| Section 5.1 | "but some it is likely the project would increase releases" | garbled, rewrite |
| Section 5.1 | "lead to reasonable validation" | "led to reasonable validation" |
| Section 4.7 | "were raised , seven of them" | "were raised, seven of them" |
| Section 5.2 | "Thirty-six ... Ten ... 5 years ... 4 where" | pick one number style |
| Section 5 | "flood risk management—it is also operated" | split into two sentences |

Four semicolons remain in body prose, in Section 4.8, the Figure 5-6 caption,
the Appendix A note and the Appendix E note.

Last, the filenames read `2026_03_09` while both documents are dated 3
September 2026. Anyone sorting by filename will read them as 9 March. Consider
`2026_09_03`.

---

## Step 12. Final verification

    cd CAS_Reg_Unreg\docs
    python "#Figure_Check.py"

Every figure should read CURRENT or BY HAND. Nothing should read DIFFERS or
NO ORIGIN.

Then tell me and I will re-run the table by table numeric check against the
CSVs, the same one that produced the review, so the final file is verified
rather than assumed.
