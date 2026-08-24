# Downstream Confluence (Coweeman + Arkansas + Ostrander) — Where This Stands

Written to consolidate a sprawling sub-task after finding that
`COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx` (uploaded to `CAS_Reg_Unreg/data/`)
already solves most of this problem. Read this before re-deriving anything
below — the precedent method is real, numeric, and closer to finished than
it looked from the CDID3 report alone.

## Script run sequence

**Group A — pre-existing pipeline, outputs already committed; only re-run to
refresh source data:**

1. HEC-SSP (external tool) -> Cowlitz unregulated Bulletin 17 report,
   `CAS_Unreg_FF/ssp/`
2. `CAS_Unreg_FF/src/Frequency_Curves_And_Table.py` -> parses #1 ->
   `CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv`
3. `CAS_Reg_Unreg/src/#Unreg_Reg_Curve.py` -> needs #2 ->
   `CAS_Reg_Unreg/output/regulated_frequency_inferred.csv`
4. HEC-ResSim "WCM_RC" run (external) ->
   `CAS_Reg_Unreg/output/ResSim_WCM_RC.dss` (both `Flow-UNREG` and `Flow`
   paths at Castle Rock)
5. `CAS_Reg_Unreg/src/#Coweeman_Timing.py` -> fetches/caches raw USGS +
   Ecology Coweeman data into `CAS_Reg_Unreg/data/coweeman/` -- **needs
   network access** (USGS hosts are blocked from this sandbox by policy;
   run locally). Cached files are already committed, so this only needs
   re-running to pull newer years.
6. `CAS_Reg_Unreg/src/#Coweeman_Proportion.py` -> needs #4 + #5's cache ->
   `coweeman_proportion.csv`

**Group B — this sub-task's new scripts, run in this order:**

7. `#Coweeman_FlowFrequency.py` -> needs #5's cached files (already
   present) -> `coweeman_frequency_table.csv`
8. `#Coweeman_RegPeak_Timing.py` -> needs #4 + #5 -- diagnostic only,
   nothing downstream depends on it
9. `#Coincident_PerfectCorrelation.py` -> needs #3 + #7 ->
   `coincident_perfect_correlation.csv`
10. `#Coincident_CorrConditioned.py` -> needs #2 + #3 + #6 + #7 (#9
    optional, only used for a comparison line on its plot) ->
    `coincident_corr_conditioned.csv`
11. `#Coincident_TieredScaling.py` -> needs #3 + #7 ->
    `coincident_tiered_scaling.csv`, the current deliverable

Since Group A's outputs are already committed, a fresh clone can run
steps 7-11 directly without touching 1-6 -- those only matter for
refreshing source data (a newer SSP run, more years of Coweeman gage
data, a different ResSim run).

## What exists today

- **Castle Rock, regulated and unregulated**: built by `#Unreg_Reg_Curve.py`
  / `CAS_Unreg_FF`'s pipeline. This curve's `reg_2009_cfs` comparison
  column, already in `regulated_frequency_inferred.csv`, turns out to be a
  literal transcription of `COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx`'s own
  Castle Rock table (Table 0 / Table B-5#) — confirmed by exact match at
  every AEP. That connection wasn't traceable to a source document until
  now; it is that document.
- **Coweeman, unregulated**: `#Coweeman_FlowFrequency.py` now adopts
  CDID3's Table 5-2 curve directly (see `CDID3_Coincident_Frequency_
  Notes.md` and the "adopt CDID3's actual curve" commit). This is a MORE
  RECENT and MORE RIGOROUS Coweeman curve than the one used in
  `COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx` (which used only 1950-1984 +
  the 1996 historic peak, no Bulletin 17C/EMA, no 2006-2015 Ecology data)
  — at 1% AEP the older document's gage-level value is 11,500 cfs vs
  CDID3's 14,100 cfs, an 18% gap. Worth carrying the better (CDID3) curve
  forward into whatever comes next, not the 2009-era one.
- **Two coincident-combination methods**, built for Castle Rock + Coweeman
  only: `#Coincident_PerfectCorrelation.py` (same-AEP sum) and
  `#Coincident_CorrConditioned.py` (correlation-conditioned, r=0.65 from
  79 concurrent events).

## The new problem

A combined regulated flow-frequency curve **downstream of Castle Rock**,
past the confluences with Arkansas Creek, Ostrander Creek, and the
Coweeman — not just Coweeman alone. Three sources of local inflow, none
of them currently gaged in a way that supports an independent Bulletin
17C curve for two of them.

## The precedent — `COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx`, Section B.6

This is the document behind "the 2009 study did X" from earlier in this
project. It solves close to exactly this problem, for these exact three
tributaries, with real numbers. Two methods are described, one superseding
the other:

**1997 method (superseded, explicitly rejected)**: drainage-area ratio
applied to the UNREGULATED discharge at Castle Rock. The document's own
words: "By using drainage area ratios, the locals below Castle Rock are
assumed to have similar rainfall intensities and identical hydrograph
shapes as the entire [Cowlitz] basin above Castle Rock. This results in an
unrealistic overestimation of locals inputs." This is essentially your
candidate **Option 1** (drainage ratios off a single curve) — already
tried, on this exact river system, and already found to overestimate.
Worth knowing before reviving it.

**2009-era method (adopted, current baseline)** — your candidate
**Option 2**, already executed:

- Arkansas Creek (44.7 sq mi) and Ostrander Creek (25.8 sq mi): no
  systematic gage record at either. Discharge-frequency statistics from
  **USGS StreamStats regional regression** (2008 version of the tool).
  Full table below.
- Coweeman: the gage-based curve described above (1950-1984 + 1996
  historic; superseded now by CDID3's curve, use CDID3's numbers instead
  of these when rebuilding this).
- **Timing**: no concurrent sub-daily data existed for these tributaries
  at the time, so timing was inferred from an ANALOG gage — East Fork
  Lewis near Heisson (USGS 14222500), a similarly-sized basin (125 sq mi)
  25 miles away. Over 17 events (WY1996-2007), the East Fork Lewis peaked
  1-12 hours before the Cowlitz's Castle Rock peak, and AT the moment of
  the Castle Rock peak, East Fork Lewis flow was 56-98% of its own peak,
  averaging **80%**.
- **Combination rule**: same-AEP pairing (like the perfect-correlation
  method), each tributary's peak scaled by the flat 80% timing factor,
  added sequentially downstream:

  ```
  below_Arkansas   = Castle_Rock_reg(AEP) + 0.80 * Arkansas(AEP)
  below_Ostrander  = below_Arkansas       + 0.80 * Ostrander(AEP)
  below_Coweeman   = below_Ostrander      + 0.80 * Coweeman(AEP)
  ```

  Verified by hand at 1% AEP: 97,000 + 0.8(4,220) = 100,376 ≈ 100,400
  (their table); +0.8(2,240) = 102,192 ≈ 102,200; +0.8(11,900) = 111,720
  ≈ 111,700. Formula confirmed exactly.

- **Its own stated caveats — read this section, it already raises your
  timing-attenuation point**: "The 80% reduction... likely overestimates
  the discharges at the smaller tributaries... especially for extremely
  large events where the regulated peak discharge at Castle Rock would
  not occur with the locals but after the reservoir fills, resulting in a
  later surge well after the other larger tributaries have peaked." That
  is exactly the mechanism you described this session — reservoir
  regulation sometimes delays the Castle Rock peak past the locals'
  peaks, sometimes doesn't, depending on whether the reservoir fills and
  starts passing inflow through before or after the local peak arrives.
  The 2009-era study flagged this as a real limitation of its own method
  and did not resolve it; it isn't a new problem you're finding, it's a
  known open one.
- Drainage-area adjustment from gage to confluence for Coweeman: 119 -> 127
  sq mi (factor 1.07). **This is the same factor CDID3 used in 2016** —
  two independent studies applying the identical correction is reasonably
  strong precedent for applying it in `#Coweeman_FlowFrequency.py` too
  (currently flagged there as an open item, unresolved). Worth closing
  that out now.
- EYR: 90 years (same as Castle Rock) applied uniformly across the whole
  reach including all three flow-change locations — no separate
  uncertainty treatment per tributary.

### Table B-VI — Arkansas Creek and Ostrander Creek (USGS StreamStats)

| AEP | Arkansas Creek (cfs), 44.7 sq mi | Ostrander Creek (cfs), 25.8 sq mi |
|---|---|---|
| 0.999 | 410 | 230 |
| 0.99 | 580 | 325 |
| 0.95 | 790 | 440 |
| 0.90 | 930 | 520 |
| 0.80 | 1,130 | 620 |
| 0.70 | 1,300 | 710 |
| 0.60 | 1,470 | 790 |
| 0.50 | 1,620 | 881 |
| 0.40 | 1,800 | 960 |
| 0.30 | 2,040 | 1,100 |
| 0.20 | 2,300 | 1,220 |
| 0.10 | 2,740 | 1,470 |
| 0.05 | 3,200 | 1,700 |
| 0.02 | 3,770 | 2,010 |
| 0.01 | 4,220 | 2,240 |
| 0.005 | 4,650 | 2,480 |
| 0.002 | 5,310 | 2,790 |
| 0.001 | 5,700 | 3,010 |

### Table B-VIII — combined result (2009-era, for reference/validation only — rebuild with CDID3's Coweeman curve before using)

| AEP | below Arkansas (RM 16.5) | below Ostrander (RM 8.54) | below Coweeman (RM 1.7) |
|---|---|---|---|
| 0.99 | 22,200 | 22,400 | 24,000 |
| 0.90 | 33,100 | 33,600 | 36,100 |
| 0.50 | 56,800 | 57,500 | 61,700 |
| 0.10 | 76,200 | 77,400 | 84,100 |
| 0.05 | 81,600 | 82,900 | 90,500 |
| 0.02 | 91,000 | 92,600 | 101,300 |
| 0.01 | 100,400 | 102,200 | 111,700 |
| 0.005 | 107,600 | 109,400 | 119,400 |
| 0.002 | 113,700 | 115,700 | 126,000 |
| 0.001 | 160,200* | 162,500* | 173,800* |

(*The table's 0.002/0.001 rows visibly jump discontinuously against the
Castle Rock curve's own jump between its 0.005 and 0.002 rows — 110,000 to
156,000 cfs, a known kink in that curve, not a transcription error here.
Carry the same caution into any rebuild.)

## How your three candidate options actually compare, given this

1. **Drainage ratios off Coweeman** — not what either precedent study did
   (1997 did ratios off Castle Rock itself and was rejected for
   overestimating; nobody tried ratios off Coweeman specifically). No
   direct precedent either way. Weakest-justified of the three unless
   there's a specific reason to prefer Coweeman over StreamStats as the
   ratio basis.
2. **StreamStats for Arkansas/Ostrander + Coweeman's own curve** — this
   IS the 2009-era precedent, already executed once. The strongest
   starting point: re-run StreamStats (the cited tool is a 2008 version;
   current StreamStats may have updated regional regression equations
   since), keep CDID3's better Coweeman curve instead of the 2009-era
   one, and decide same-AEP-plus-80% vs. something more like the
   correlation-conditioned approach already built for Coweeman.
3. **A single flat scaling factor for all three** — functionally what the
   2009-era 80% already is, just derived per-basin-implicitly rather than
   applied as one number across all three. If the goal is literally "one
   scale factor, applied to all three," the 80% figure is sitting right
   there, empirically derived (East Fork Lewis analog), real precedent —
   though derived from a proxy basin, not these basins directly.

**Given the precedent, Option 2 is the natural default** — not because
the other two are wrong, but because it's the one with a worked example,
real StreamStats numbers, and an explicit, already-catalogued list of its
own weaknesses to either accept or improve on. That's a decision for you
to make, not something to default into silently.

## Two things worth fixing before adopting any of this outright

### 1. The correlation estimate (r=0.65) may itself be contaminated by rating-curve capping

Checked directly: of the 79 concurrent events behind
`coweeman_proportion.csv` and `#Coincident_CorrConditioned.py`'s r=0.65,
**32 of 79 (41%) have their Coweeman peak reading flagged quality code 10
("above rating, but within 2x")** — the Ecology gage's rating curve is
being extrapolated, not measuring directly. Most concerning: **the single
largest event in the whole 79-event record — Dec 9, 2015, 152,270 cfs at
Castle Rock — is one of them.** Its reported Coweeman value (3,328 cfs) is
almost certainly an understatement; this is the identical Dec 2015 event
CDID3's own report says has "no reliable flow measurement" at all. Of the
9 events exceeding 60,000 cfs at Castle Rock (the bin that matters most
for the flood-tail correlation), 2 carry this flag — a smaller fraction,
but one of the two is the single most influential point in the sample.

This doesn't necessarily mean r=0.65 is wrong, but it means the number is
built partly on capped values at exactly the tail where it matters most,
and the direction of the bias isn't obvious without more work (capping
compresses the top of the Coweeman range, which could either weaken the
apparent correlation or just add noise, depending on how consistently
different storms hit the same ceiling). Worth deciding whether to accept
this as another "interim, not final" caveat (consistent with everything
else in this analysis) or to try correcting the worst-affected events
using an analog-gage approach similar to what the 2009-era study did for
Arkansas/Ostrander timing — you already have real concurrent data now
that didn't exist in 2009, which is a genuine opportunity to do this
better than either precedent could.

### 2. Coweeman gage-to-confluence scaling (119 -> 127 sq mi, factor 1.07)

Flagged as an open item in `#Coweeman_FlowFrequency.py`. Now precedented
twice (CDID3 2016, this document 2009-era) — reasonable to just apply it
rather than continue flagging it as unresolved.

### 3. The 2009-era method pairs a REGULATED Castle Rock AEP against NATURAL tributary AEPs — a real, separate conflation

Caught when asked directly whether the paired Castle Rock value was
regulated or unregulated. It's regulated (Section B.6 title, and the
Castle Rock column is an exact match to `reg_2009_cfs`). Whether that's a
coherent thing to pair against three naturally-driven tributary curves
depends on how the regulated curve's AEP axis was built, and there are two
different cases:

- **Transform-based** (rank-preserving by construction): each unregulated
  annual value is routed through the dam and KEEPS its original rank, so
  regulated-AEP and unregulated-AEP are literally the same axis, read off
  two different y-curves. This is how **our own** Castle Rock pipeline
  works (`#Unreg_Reg_Curve.py`: `reg_inferred_cfs` is a transform of
  `unreg_computed_cfs` at the same AEP, one shared AEP column). Under this
  construction, pairing tributary AEP against regulated or unregulated
  Castle Rock AEP is identical -- no conflation.
- **Independently, graphically fit** to the regulated series' own rank
  order: this is what `COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx` did --
  "the regulated discharge-frequency curve at Castle Rock is graphically
  fit," "graphically drawn... through the 83 years" (paras 14, 48) --
  fit directly to the regulated annual-peak series, not derived by
  transforming the unregulated curve. Regulation can and does reorder
  years (a middling natural year can produce a relatively large regulated
  peak if the reservoir was already full going in; we've seen this same
  rank-sensitivity in the reg/unreg transform work elsewhere in this
  project, e.g. the transform slope exceeding 1 near pass-through). Under
  this construction, "Castle Rock regulated 1% AEP" and "Arkansas Creek
  natural 1% AEP" are not guaranteed to represent the same underlying
  storm year.

This document is the second case. So pairing its regulated Castle Rock
column against Arkansas/Ostrander/Coweeman's natural AEPs is a real,
separate weakness from the flat-80%-timing-factor caveat already
catalogued above -- worth listing on its own.

**Practical upshot**: this doesn't apply to a rebuild on our own current
Castle Rock curve, which IS transform-based (case 1). One more reason,
on top of the better Coweeman curve and current StreamStats numbers, to
rebuild the downstream-confluence combination fresh rather than splice
the 2009-era regulated-Castle-Rock column into a new combination.

### 4. Is the borrowed 80% timing factor actually evidenced? Split answer -- no for Arkansas/Ostrander, and for Coweeman the real number is lower at the tail

Arkansas Creek and Ostrander Creek: no, and it can't be -- the document
says so directly ("no long-term, systematic records of peak flow available
for Arkansas Creek and Ostrander Creek"), which is exactly why it reached
for the East Fork Lewis analog instead of measuring either creek directly.
Nothing in this repo changes that; worth a check for a newer gage but not
assumed.

Coweeman: unlike in 2009, we now HAVE the concurrent data to check this
directly, instead of borrowing an analog basin. Computed straight from
`coweeman_proportion.csv` (79 events): the fraction of the Coweeman's OWN
peak it is running at when the (unregulated) Castle Rock peak arrives --
the same quantity the East Fork Lewis 80% describes, for the Coweeman
itself.

| | n | mean | median |
|---|---|---|---|
| Overall | 79 | 0.740 | 0.775 |
| 20,000-40,000 cfs bin | 51 | 0.794 | 0.809 |
| 40,000-60,000 cfs bin | 19 | 0.732 | 0.764 |
| **>60,000 cfs bin** | **9** | **0.457** | **0.515** |
| 2009-era analog (East Fork Lewis, n=17) | -- | 0.80 | -- |

At moderate magnitude the borrowed 80% actually holds up well -- the
20-40k bin's median (0.809) nearly matches it. At the flood tail, where a
frequency study lives, the real Coweeman-specific number is roughly HALF
the borrowed one (0.457-0.515, not 0.80). This is the same mechanism the
2009-era report's own caveat gestures at conceptually ("especially for
extremely large events... a later surge well after the other larger
tributaries have peaked") -- now with a number behind it, for one of the
three creeks.

**Checked the overlap directly rather than leaving it as a caveat** (asked
directly: how do you know the tail isn't just the capping?). Pulled the
quality code on BOTH readings behind each of the 9 tail-bin ratios (the
"own peak" denominator and the "value when Castle Rock peaks" numerator):

| event_time | Castle Rock (cfs) | own peak (code) | at-Castle-Rock-peak (code) | ratio |
|---|---|---|---|---|
| 2006-11-07 | 155,018 | 7,297.5 (160) | 4,062.5 (160) | 0.557 |
| 2015-12-09 | 152,270 | 3,327.5 (10) | **254.0 (MISSING)** | 0.076 |
| 2009-01-08 | 143,468 | 7,810.0 (160) | 2,442.5 (160) | 0.313 |
| 2007-12-04 | 89,749 | 5,292.5 (160) | 2,780.0 (160) | 0.525 |
| 2011-01-17 | 85,794 | 5,372.5 (100) | 2,765.0 (10) | 0.515 |
| 2015-11-18 | 73,658 | 6,120.0 (100) | 2,422.5 (10) | 0.396 |
| 2017-03-16 | 68,453 | 3,335.0 (10) | 2,556.0 (**MISSING**) | 0.766 |
| 2006-12-15 | 63,582 | 3,025.0 (160) | 2,127.5 (50) | 0.703 |
| 2008-11-13 | 62,213 | 7,640.0 (160) | 2,005.0 (10) | 0.262 |

One event (2015-12-09) is not a capped-but-real reading -- it's
interpolated across an actual data gap (`MISSING`), during the exact
event CDID3's own report says has no reliable flow measurement at all.
That one is fair to drop outright, no ambiguity. Its ratio (0.076) is by
far the lowest in the bin, so it looked like it might be doing a lot of
work.

It isn't, though: dropping it moves the bin median from 0.515 to 0.520 --
essentially nothing. Also dropping the other `MISSING`-flagged point
(2017-03-16) leaves median 0.515 on n=7. The pattern -- tail running well
below 0.80 -- does NOT depend on the one point that was genuinely broken.

What it DOES depend on, and what I should be more careful claiming: the
remaining seven readings are not all clean. Their codes are 10 ("above
rating, within 2x" -- a real reading, modestly extrapolated), 50
("estimated"), 100 ("modeled flow"), and an undocumented code **160** that
shows up on several of the oldest events and that I have not found a
legend definition for anywhere in the cached files. I don't know how much
error 160 represents. So the correct claim is narrower than "the tail
ratio is 0.46-0.52": it's "the tail ratio is well below 0.80, and that
conclusion survives removing the one point known to be invalid, but most
of the remaining points still carry non-good quality flags of uncertain
size." Real and worth acting on, not fully resolved -- same posture as
everything else flagged in this document, not a clean result.

Also worth remembering: this is measured against UNREGULATED Castle Rock
timing (matching how `coweeman_proportion.csv` was built), not the real
regulated gage timing the East Fork Lewis figure used -- a true
apples-to-apples check would use the regulated peak timestamps in
`coweeman_event_timing.csv` instead, not yet done.

**Implication if this holds up**: using the flat 80% for the Coweeman at
flood-tail AEPs would overstate its coincident contribution to the
combined downstream peak by roughly a factor of 1.5-1.7x versus what the
Coweeman's own concurrent record actually shows. (The mechanism proposed
here originally -- reservoir regulation delaying the Castle Rock peak
past the locals' peaks -- turns out NOT to be well supported once tested
directly; see #5 below. The tail drop-off itself is real, but it isn't
because the regulated peak arrives meaningfully later than the
unregulated one.)

### 5. Everything above uses UNREGULATED Castle Rock timing -- tested against the REGULATED peak directly, and the "regulation delays it" story doesn't hold up

Asked directly: all of #4 measures the Coweeman against the moment the
*unregulated* Castle Rock peak arrives, not the *regulated* one -- and the
regulated peak is the one that actually matters for a downstream combined
flow. If regulation delays the peak (the mechanism assumed above), the
Coweeman should have receded further still by the time the real,
regulated peak shows up.

Tested this directly with a new script, `#Coweeman_RegPeak_Timing.py`,
using the SAME ResSim run's regulated companion series
(`//CastleRock_NWS/Flow//1Hour/ResSim_WCM_RC/`, current operating rules
applied to the same reconstructed unregulated inflow -- the right,
consistent basis, matching how the rest of this project's regulated
curves are built, though NOT the actual historically-observed regulated
record for these specific events).

**Lag (regulated peak time minus unregulated peak time)**:

| bin | n | mean | median | range |
|---|---|---|---|---|
| 20-40k | 51 | -0.6h | +0.0h | -72h to +37h |
| 40-60k | 18 | -0.5h | +0.0h | -13h to +7h |
| **>60k** | **8** | **-3.2h** | **-3.5h** | **-15h to +5h** |

The lag is small at every magnitude, and at the flood tail it's slightly
NEGATIVE on average -- the regulated peak tends to arrive a few hours
EARLIER than the unregulated one in this ResSim run, not later. That's
the opposite of the "reservoir fills and passes flow through late"
mechanism assumed above (and stated as a caveat in the 2009-era report
itself). It can still happen event-by-event -- the range reaches +37h in
the smallest bin -- but it isn't the dominant pattern in this data.

**Ratio, unregulated-timing vs. regulated-timing**:

| bin | n | ratio @ unreg-timing (mean/median) | ratio @ reg-timing (mean/median) |
|---|---|---|---|
| 20-40k | 51 | 0.792 / 0.809 | 0.791 / 0.809 |
| 40-60k | 18 | 0.721 / 0.762 | 0.730 / 0.762 |
| **>60k** | **8** | **0.533 / 0.520** | **0.571 / 0.494** |

Switching to regulated-peak timing barely moves the numbers -- mean goes
up slightly, median down slightly, both still well below 0.80. **The
tail-drop-off finding in #4 survives this test**, but the explanation for
WHY changes: it is not primarily a timing-lag effect. Something else is
driving it (possibly the largest storms' hydrograph shape putting the
Coweeman's own peak and Castle Rock's peak on genuinely different points
of two different-shaped curves, independent of any regulation-driven
delay) -- not yet investigated further.

**Two things caught while building this, flagged rather than smoothed
over**:
- One event (2017-03-16) shows ratio=0.9965 here, at both timings --
  but the direct quality-code check in #4 found this same event's
  at-Castle-Rock-peak reading flagged `MISSING` in the raw file. The two
  scripts are reading the same underlying data two different ways
  (exact-timestamp lookup in #4 vs. `.resample('1h').mean()` here), and
  the resample appears to be interpolating across the same gap the exact
  lookup correctly flagged as missing. This is the same underlying
  problem wearing a different disguise -- worth fixing how gaps are
  handled before trusting this event either way.
- 1 of 77 events had its regulated peak land within 3 hours of the
  +/-72 hour search window edge -- flagged by the script itself, not
  yet re-run with a wider window to confirm the true peak was found.

## Not yet decided — flagging rather than picking

- Whether to rebuild the downstream-confluence combination fresh (own
  StreamStats pull, own timing analysis) or start from the 2009-era
  numbers and only replace the Coweeman piece.
- Whether "timing" for Arkansas/Ostrander should stay as the 2009-era
  flat 80% (proxy-basin-derived) or get its own analog-gage or
  correlation treatment the way Coweeman now has.
- Whether the capped-event correlation issue needs fixing now or can
  ride as a documented caveat, same posture as the interim Coweeman
  curve before CDID3's was adopted.
