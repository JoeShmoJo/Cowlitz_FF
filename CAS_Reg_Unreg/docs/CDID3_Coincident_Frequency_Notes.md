# CDID3 Phase 2a Report — Coincident Frequency Notes

Compacted from `CAS_Reg_Unreg/data/CDID3_Phase2a_Report_Final_2016-06-30.docx`
(USACE Portland District, final 30 Jun 2016) so the full 8.4 MB source doesn't
have to stay in context. Read the source docx again only if a number or a
methodological detail below turns out not to be enough — it's still in the
repo.

## What the report is

A FEMA-accreditation levee study (EC 1110-2-6067 / 44 CFR 65.10) for the
CDID3 levee along the **lower Coweeman River** near its confluence with the
Cowlitz at Kelso. The Coweeman levee segment's peak water surface elevation
is driven by **three partially-correlated flooding sources**: Coweeman flow
(unregulated), Cowlitz flow (regulated, Mossyrock/Mayfield since 1968), and
Columbia River stage (tidal, at Longview). This is the same structural
problem our study has — an unregulated tributary's peak, a regulated
mainstem peak, and a question about whether/how to combine them — except
CDID3 solves it downstream of Castle Rock and includes a third (tidal)
source we don't have.

**Scope limit stated explicitly in the report**: results are for the 1% ACE
event only. Bulletin 17C. No general stage-frequency curve is produced,
and the authors say not to extrapolate to more frequent events (curve
would be wrong below the 1% ACE tail — the seasonal-max approximation
breaks down there, see below).

## Why winter-only is already their assumption, not just ours

Section: *Seasonality*. They collapse the annual analysis to **winter
(Nov–Mar) only**, same as we want, and give a concrete justification we can
reuse in a memo: at Castle Rock 1927–2014, 82 of 88 annual peaks fell
Nov–Mar; spring/summer Coweeman and Cowlitz frequency curves are far lower;
and even though the Columbia has a spring freshet peak, the *winter* Columbia
curve is close to the annual curve for large events while the *spring* curve
is much lower. Net effect: for large events the annual result is dominated by
winter, so treating winter as the surrogate for annual loses essentially
nothing at the tail — but they flag (as we should) that this substitution is
NOT valid for common events (e.g. 50% ACE), only for the flood tail.

They also state the implicit assumption plainly: using a winter-season-max
frequency curve for each source assumes the winter maximum of all three
sources tends to come from the *same storm* — supported in their record by
named concurrent events (Feb 1996, Jan 2009, Dec 2015).

## Their regulated Cowlitz curve — already exactly what we need

Section: *Cowlitz River Flow Frequency Curve*. They did NOT derive a new
Cowlitz curve — they took the **regulated** flow-frequency curve straight
from Table 3-1 of the 2009 Level of Protection report (`USACE 2010b`),
treated as representing the winter season (>90% of annual peaks fall
Nov–Mar), with **90 equivalent years of record** assumed for its confidence
limits. Source blend behind that 2009 curve: observed 1969–2009 (post-dam),
simulated regulated flows for large pre-1969 events, an unregulated→regulated
conversion for other pre-1969 peaks, and synthetic hydrology peaks. This is
the same lineage `CAS_Reg_Unreg/src/#Unreg_Reg_Curve.py` is extending/
reworking now — worth citing as precedent that literal use of a regulated
curve (not re-derived from an unregulated one) is an accepted approach here,
and that 90 EYR is the number USACE has already put in front of a reviewer
for this exact regulated curve at a nearby location.

## Their Coweeman (unregulated) curve

Section: *Coweeman River Flow Frequency Curve*. Bulletin 17C draft methods
(pre-publication in 2016), Log-Pearson III, EMA with multiple Grubbs-Beck,
regional skew blended in (0.2, MSE 0.112, Cooper 2005). Perception
thresholds handled historic/censored data explicitly — worth reusing as a
template if we ever formalize censored years:

| Start | End | Low threshold (cfs) | High threshold |
|---|---|---|---|
| 1933 | 1949 | 12,000 | ∞ |
| 1950 | 1984 | 0 | ∞ |
| 1985 | 2006 | 11,700 | ∞ |
| 2007 | 2015 | 0 | ∞ |
| 2016 | 2016 | 11,700 | ∞ |

44 equivalent years of record assumed. Drainage-area-ratio scaled from gage
(119 mi²) to mouth (127 mi²), factor 1.07. They explicitly *rejected* a
two-station comparison with East Fork Lewis (best-correlated neighbor gage)
because concurrent vs. non-concurrent mean/stdev differed <5% — not worth the
complexity/uncertainty it would add. (Analogous to our own
`#Coweeman_Proportion.py` reasoning about not injecting extra uncertainty
through indirect transforms.)

## Full frequency-curve table (Table 5-2, winter season, all three sources)

| AEP | Coweeman unreg (cfs) | Cowlitz reg (cfs) | Columbia stage (ft NAVD88) |
|---|---|---|---|
| 99.99% | 1,700 | 10,400 | 9.0 |
| 99% | 2,500 | 11,700 | 9.6 |
| 95% | 3,100 | 19,000 | 10.8 |
| 90% | 3,500 | 25,200 | 11.4 |
| 80% | 4,000 | 33,400 | 12.2 |
| 70% | 4,500 | 38,100 | 12.7 |
| 60% | 4,900 | 42,800 | 13.2 |
| 50% | 5,400 | 48,000 | 13.6 |
| 40% | 6,000 | 53,200 | 14.4 |
| 30% | 6,600 | 60,500 | 15.3 |
| 20% | 7,500 | 68,800 | 16.4 |
| 10% | 9,000 | 83,400 | 17.8 |
| 5% | 10,500 | 99,900 | 18.8 |
| 2% | 12,500 | 112,600 | 20.5 |
| 1% | 14,100 | 118,200 | 21.3 |
| 0.5% | 15,800 | 129,700 | 22.7 |
| 0.2% | 18,200 | 166,500 | 24.3 |
| 0.1% | 20,100 | 197,000 | 25.0 |
| 0.05% | 22,200 | 307,500 | 26.1 |
| 0.01% | 27,400 | 398,800 | 28.4 |

(Cowlitz column is the 2009 LOP regulated curve reproduced verbatim; note its
1%-AEP value of 118,200 cfs sits well below our `BAND_CLIP_MIN_CFS = 60,000`
physical-ceiling discussion territory only at the low end — at the 0.5–0.01%
tail it's clearly pass-through-dominated, consistent with what we've already
found about the transform's slope exceeding 1 near storage exhaustion.)

## Coincidence methodology — this is the part that answers your question

Section 6, *Coincidence Analysis Inputs*. They explicitly reject both
extremes (fully correlated = old, over-conservative Phase 1 approach; fully
independent = unrealistic) and instead separate the joint-event problem into
two independent axes, both handled by Monte Carlo, not by a single offset
correction:

### 1. Event magnitude — correlation, not a fixed ratio

Historic concurrent events (peaks within a few days of each other) are each
converted to a z-statistic via their own frequency curve (i.e., the ACE of
that event on that source's curve, normal-transformed), and a Pearson
correlation matrix is built across the three sources from those z's:

| | Coweeman | Cowlitz | Columbia |
|---|---|---|---|
| Coweeman | 1 | 0.814 | 0.327 |
| Cowlitz | 0.814 | 1 | 0.358 |
| Columbia | 0.327 | 0.358 | 1 |

(Final adopted matrix, after adjusting the Coweeman-Columbia cell down to be
consistent with the Cowlitz-Columbia full-period-of-record correlation,
since Coweeman-Cowlitz stays essentially flat between concurrent-period and
full-period estimates — see docx §6.1 for the adjustment formula if it's
ever needed.) **Cowlitz-Coweeman correlation is strong (0.81) — adjacent
watersheds, same storms — but not 1.0**, because storm type shifts the ratio
(a short intense cell can spike the small Coweeman basin harder in ACE terms
than the large regulated Cowlitz). That imperfect correlation is exactly the
"offset" your question is asking how to acknowledge: their answer is not a
static offset but a sampled joint distribution.

### 2. Event timing — representative shape sets, not a fixed lag

Separately from magnitude, they characterize *when* each source peaks
relative to the others using a box-and-whisker of observed lag (Coweeman
first, Cowlitz next, Columbia last, on average — but "not always," with real
spread). Rather than model timing continuously, they pick **three
representative historic events as timing templates**, sampled with unequal
weight in the Monte Carlo, and never mixed (a whole event's timing shape is
used together):

| Scenario | Event | Weight | Character |
|---|---|---|---|
| A | Nov 2012 | 3 | Cowlitz/Coweeman peak within 4 hrs (near-coincident) |
| B | Feb 2012 | 2 | Cowlitz/Coweeman peak within 8 hrs |
| C | Jan 2011 | 1 | Columbia peaks much later |

### 3. Monte Carlo combination (HEC-WAT)

Two-loop sampling: outer loop = knowledge uncertainty (sampled once per
*realization*, e.g. which point on the frequency-curve confidence band),
inner loop = natural variability (sampled once per *event*: magnitude via
the correlation matrix + a timing shape set). 100 realizations × 500 events
(50,000 total) fed through a calibrated HEC-RAS model gives a simulated peak
WSE distribution at the levee; the 99th percentile of each realization is
that realization's 1% ACE estimate, and the spread of those 100 estimates
across realizations gives the hydrologic uncertainty (checked and accepted
as normally distributed via K-S test, p=0.25). Convergence was checked by
percent-error-of-the-mean formula, not by eyeballing a single stabilization
plot — worth citing if a reviewer asks how we know our own Monte Carlo
(`RUN_MONTE_CARLO_CHECK` in `#Unreg_Reg_Curve.py`) is converged.

### Total uncertainty — RSS combine, same shape as our own

Section 8. Hydrologic uncertainty (from the coincidence Monte Carlo, ranges
1.6–2.4 ft across the levee, worse near the confluence than upstream — makes
sense, that's where all three sources' relative effects are strongest) is
combined via **root-sum-of-squares** with hydraulic natural uncertainty
(Manning's n sensitivity, treated as knowledge uncertainty since the model
can't sample it directly — this is the same "n can't be sampled directly in
Monte Carlo, treat as a fixed knowledge term instead" workaround, same
justification we've been using), hydraulic model uncertainty (0.7 ft, from
EM 1110-2-1619 Table 5-2 given fair calibration confidence), and sediment
uncertainty (0.7 ft flat, borrowed directly from the 2009 LOP report given
ongoing Mount St. Helens/Toutle sediment supply effects on the Cowlitz
channel — this reach is literally just downstream of ours, so 0.7 ft is a
directly reusable number if we ever need a Castle-Rock-adjacent sediment
term). This is the same EM 1110-2-1619 Sec. 4-4.b(3)/Eq. 4-6 RSS combination
already being used in `#Unreg_Reg_Curve.py` — same standard, applied one
level up (stage, not flow, and combining a Monte-Carlo-derived hydrologic
term with fixed hydraulic terms rather than combining two flow-space terms).

## The open question this doesn't answer for us

CDID3 solves the coincidence problem **downstream** of Castle Rock, where the
quantity of interest is a *stage* produced by routing all three inputs
through a shared hydraulic model (HEC-RAS) — the "offset" is absorbed by
that model run, not by an analytic adjustment to any one curve. Our problem
is upstream of that: we want a coincident **regulated peak flow-frequency**
at (or near) Castle Rock itself, driven by Cowlitz regulated flow and
Coweeman unregulated flow, without necessarily running every sampled pair
through a hydraulic/routing model the way they did. Two live options, not
yet decided:

1. **Full replication of their approach** at our own confluence: build a
   correlation matrix (Cowlitz-Coweeman only, drop Columbia — not relevant
   this far upstream) from concurrent events using `#Coweeman_Timing.py` /
   `#Coweeman_Proportion.py`'s existing event-matching, pick representative
   timing shape sets from the same two scripts' output, and Monte Carlo
   sample joint (magnitude, timing) pairs the way HEC-WAT does — but without
   a hydraulic model to route them through, since we want a flow answer, not
   a stage answer. That simplification (skip the routing step) is actually
   *easier* than what CDID3 did, since routing was their reason for needing
   HEC-RAS at all.
2. **A lighter acknowledgment-only treatment**: use the existing
   `#Coweeman_Proportion.py` ratio_coincident / ratio_peak split (already
   quantifies exactly the "timing penalty" concept CDID3's shape-set weights
   encode) as a documented, memo-level caveat on the combined peak rather
   than folding it into the frequency curve numerically — i.e., report the
   regulated Castle Rock curve as-is and separately bound how much Coweeman
   could add if coincident, sourced from the already-built proportion/timing
   analysis, without building a new Monte Carlo layer.

`#Coweeman_Timing.py` and `#Coweeman_Proportion.py` already do the empirical
groundwork CDID3's Section 6.1/6.2 does (event magnitude ratio, event timing
lag) using our own gages, independently of this report — the CDID3 method is
useful mainly as the **precedent for how USACE has formally combined this
kind of ratio-plus-lag information into a defensible frequency product**,
not as a source of numbers we're missing. Worth deciding which of the two
options above before writing any code.
