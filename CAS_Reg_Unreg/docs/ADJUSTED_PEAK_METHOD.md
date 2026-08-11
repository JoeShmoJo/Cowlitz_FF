# Adjusted Castle Rock Peak Record — Method

Produced by `src/#Adjusted_Peak_Record.py`. Written 11 Aug 2026.

## What this does

Two ResSim runs share the same observed hydrology and the same current
operating rules, and differ only in the reservoir's starting pool:

| Run | Starting pool | Coverage |
|---|---|---|
| `ResSim_WCM_RC` | WCM rule curve | 98 water years, Oct–May windows |
| `ResSim_Obs_RC` | Observed elevation at event onset | 53 water years, 31-day windows |

Because operations are held constant, the difference between them at Castle
Rock isolates the effect of the starting pool alone. Where `WCM_RC > Obs_RC`,
starting at the rule curve produced a **higher** regulated peak than starting
from the observed pool — meaning the observed pool held storage the rule curve
does not credit, and the historical record benefited from it.

The adjustment removes that benefit from the observed record:

```
adjusted_peak = usgs_peak + (wcm_peak − obs_peak)     when wcm > obs
adjusted_peak = usgs_peak                             otherwise
```

**The adjustment is one-sided by design: the observed peak is only ever
increased, never reduced.** A negative difference means the observed pool
started *above* the rule curve, so the historical operation was already at
least as good as following the WCM — there is nothing to remove and the
observed peak stands. Only WY1976 falls in this category (Obs_RC 103,012 vs
WCM_RC 75,597; the pool started 7.5 ft above the rule curve).

The result is an observed-magnitude peak record on a consistent rule-curve
starting-pool basis.

## Peaks are matched by event, not by water year

`PEAK_METHOD = "event"` takes each run's peak **within ±3 days of the observed
peak date**, rather than each run's water-year maximum.

This matters. The Obs_RC run simulates a 31-day window around each year's
WCM_RC peak. A low starting pool attenuates the first storm more, so a *later*
storm in the same window can become that run's annual maximum. Differencing
annual maxima would then subtract two different storms.

The effect is large: event matching yields 40 usable years, annual maxima
yield 29.

## Screening

A year must pass both screens to be adjusted:

1. **Timing** — all three peak times within `EVENT_WINDOW_DAYS` (3) of each
   other.
2. **Containment** — simulated data must exist within that window. If the
   Obs_RC 31-day window does not cover the observed peak date, the runs never
   simulated that storm.

Years failing either screen are **carried through unadjusted**, not dropped, so
the record stays complete and every decision is recorded in `decision` in the
output CSV. Set `DROP_FAILED_SCREEN = True` to omit them instead.

## Results (51 shared water years, 1974–2024)

| Outcome | Years |
|---|---|
| Adjusted (increased) | 39 |
| Not adjusted — no gain from the rule curve | 5 |
| Not adjusted — failed screening | 7 |

Adjustment where applied: median **+10,966 cfs**, range +263 to +24,013 —
a median **+21.6%** on the observed peak. No year is reduced.

### The seven screened out

All seven fail for the same reason: the Obs_RC window does not cover the
observed peak date.

| WY | Observed peak | Note |
|---|---|---|
| 1977 | 15,000 | low peak |
| 1980 | 97,000 | **Mount St. Helens** |
| 1984 | 49,600 | |
| 1994 | 17,300 | low peak |
| 1995 | 47,500 | |
| 2001 | 11,600 | low peak |
| 2019 | 26,000 | |

**WY1980** is also in `MANUAL_EXCLUSIONS`. The observed 97,000 cfs peak on
18 May 1980 is the Mount St. Helens eruption — a debris-flow event with no
meteorological analogue. Both simulations peak on 18 Dec 1979 at 38,587 cfs
from an ordinary winter storm, 152 days earlier. It is excluded explicitly so
the reason is recorded rather than inferred from a date gap. Whether WY1980
belongs in a flood frequency analysis at all is a separate question this
script does not decide.

Three of the remaining six are low-peak years (< 20,000 cfs), where the annual
maximum is often not a distinct storm. WY1984, WY1995 and WY2019 are moderate
peaks worth reviewing individually — the Obs_RC window for those years was
placed on a different storm than the one that produced the observed peak.

## Two things to check before using this record

**A cluster at +14,000 cfs.** 13 of the 39 adjustments land within ±50 cfs of
+14,000, and more within ±1,000. This is most likely a fixed release rule (a
channel capacity or maximum release constraint) binding in one run and not the
other, rather than a pure starting-pool effect.

Reviewed and accepted (Aug 2026): an artifact of this size is not significant
against a 14,000 cfs adjustment, and the record is used for a
regulated-vs-unregulated relationship at large events, not for an analytical
frequency fit. The script still prints the warning so the cluster stays visible.

**Seven adjustments exceed 40% of the observed peak**: WY1989 (+59%),
WY1992 (+52%), WY1993 (+95%), WY2005 (+47%), WY2010 (+52%), WY2014 (+43%),
WY2023 (+41%). All are moderate-peak years where a large absolute adjustment
is a large fraction. WY1993 nearly doubles a 14,700 cfs peak. These are
mechanically correct given the inputs, but they move the low end of the
frequency curve substantially and deserve a look.

**WY1976** is the one negative difference: Obs_RC peaks at 103,012 cfs against
WCM_RC's 75,597, because the pool started 7.5 ft *above* the rule curve. No
adjustment is applied — the record is increase-only.

## Outputs

| File | Contents |
|---|---|
| `output/adjusted_peaks.csv` | every shared year: three peaks, delta, adjusted peak, screening result, `decision` text |
| `output/adjusted_peaks_ssp.csv` | `WY, adjusted_peak` for HEC-SSP import |
| `output/adjusted_peaks.dss` | `/COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW-ANNUAL PEAK-ADJUSTED//1DAY/ADJUSTED/` |
| `output/diagnostics/adjusted_peaks.png` | three records overlaid, plus the adjustment applied |
| `output/diagnostics/adjusted_peaks_screening.png` | timing spread per year against the screen |

The DSS record is written as a regular **1DAY** series, not IR-CENTURY: handing
an IR-CENTURY pathname to `put_ts` with a regular container **segfaults** this
pydsstools build rather than raising. The peak sits on its observed peak date;
every other day is the missing sentinel, so it plots as points in DSSVue.

## Intended use

This record feeds a **critical duration analysis** relating unregulated AEP
flows to regulated flows at large events. It is deliberately NOT intended for
an analytical frequency fit on regulated flows, which do not follow an
analytical distribution. That is why the screening losses are acceptable: the
seven excluded years are dominated by low peaks that do not influence the
large-event relationship.


## Settings that change the answer

| Setting | Default | Effect |
|---|---|---|
| `PEAK_METHOD` | `"event"` | `"annual"` uses water-year maxima instead — 29 usable years |
| `EVENT_WINDOW_DAYS` | 3 | Both the matching half-width and the timing screen |
| `REQUIRE_IN_OBS_WINDOW` | `True` | Containment screen |
| `DROP_FAILED_SCREEN` | `False` | `True` omits failed years instead of carrying them unadjusted |
| `MANUAL_EXCLUSIONS` | `{1980: ...}` | Years never adjusted, with the reason recorded |
