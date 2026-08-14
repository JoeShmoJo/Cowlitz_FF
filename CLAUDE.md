# Notes for Claude working in this repo

Written by Claude, for Claude. Everything here was verified by running it in this
repo, not assumed. Read this before telling the user what you can or cannot do.


## DSS: you CAN read and write version 6 on Linux

This is the single most important item. An earlier session got it wrong and told
the user, incorrectly, that DSS-6 files were unreadable. They are not.

    HecDss.Open(path, version=6)   # works -- reads AND writes v6 on Linux
    HecDss.Open(path)              # FAILS on a v6 file

The `DSS version 6 is not supported on Unix(Linux/MacOS)` error, and the
`DssStatusException (-123, 'DSS6 not supported in Mac and Linux')` it raises,
only fire during **auto-detection** -- that is, when the `version` argument is
omitted and pydsstools inspects the header itself. Passing `version=6`
explicitly bypasses the guard entirely. The DSS-6 Fortran layer (`zopen6Int.f`,
`ZOPEN6`) is compiled into the Linux `.so`; it is a refusal, not an absence.

Verified: created a new v6 file and wrote a record to it; opened the 55 MB
`CAS_Unreg_FF/data/obsData.dss` and listed all 9,686 records.

If a DSS call fails, try the explicit `version=` argument before concluding
anything about capability. Do not tell the user a file is unreadable until you
have tried both 6 and 7.


## Which files here are which version

    v6   CAS_Unreg_FF/data/obsData.dss          (the shared observed store)
    v6   CAS_Unreg_FF/output/MOS_Cleaned.dss
    v6   CAS_Unreg_FF/output/CAS_Unreg_SSP.dss
    v6   CAS_Unreg_FF/ssp/2026_Restudy/2026_Restudy.dss
    v6   CAS_Reg_Unreg/output/ResSimInflows.dss
    v7   CAS_Unreg_FF/ssp/2026_Restudy/Bulletin17Results/*/*.dss
    v7   CAS_Reg_Unreg/output/diagnostics/MOS_Special_Release.dss

Check the header directly when unsure -- byte 12 of the file is the version:

    open(path,'rb').read(16)   # b'ZDSS...\x06\x00\x00\x00' = v6


## pydsstools gotchas that have actually bitten

- `put_ts` needs a positional `TimeSeriesContainer`. The keyword form
  (`d.put_ts(path, values=..., interval=...)`) raises `TypeError`.
      tsc = TimeSeriesContainer(pathname, count, interval, values=...,
                                start_time="01JAN2025 0100",
                                data_units="CFS", data_type="INST-VAL")
      d.put_ts(tsc, store_flag=2)
- `ts.times` yields `HecTime` objects in some builds and **plain strings** in
  others (the user's Windows build returns strings; this sandbox returns
  HecTime). `h.datetime()` therefore raises `AttributeError: 'str' object has
  no attribute 'datetime'` on their machine. `ts.startDateTime` is `None` here
  but populated there. Build the index from the FIRST stamp plus `ts.interval`
  (seconds) instead of walking every element -- it is version-safe and far
  faster on a 859k-value record. See `first_stamp` / `series_step` in
  `#Create_Unreg_Ensembles.py`.
- `ts.times` is a generator: `len()` fails, so use `next(iter(ts.times))`.
- The CATALOG API differs between builds too. This sandbox has `search_path`;
  the user's build does NOT and raises
  `AttributeError: 'Open' object has no attribute 'search_path'`. Their build
  has `path_dict(pattern)`, which returns a dict of record-type -> list of
  pathnames (the regular time series are under `'ts-reg'`). Use the
  `catalog_paths()` helper in `#Extract_Ensemble_To_Timeseries.py`, which tries
  `search_path`, then `path_dict`, then `getPathnameList`. Never call a
  pydsstools catalog method directly.
- DSS uses **midnight-as-2400**, so a string stamp can be `01Oct1973 24:00:00`,
  which means 02Oct1973 00:00. `pd.Timestamp` raises
  `DateParseError: hour must be in 0..23` on it. Rewrite the hour to 00 and add
  a day -- see `first_stamp` in the ensemble scripts. Daily records are the ones
  that hit this.
- `ts.numberValues` is `None` here; use `len(ts.values)`.
- Missing values: mask with `np.array(ts.nodata, dtype=bool)`, and also treat
  `<= -900` as missing. Write missing as `-3.4028234663852886e38`.
- `search_path('/*/*/*/*/*/*/')` is the catalog call. There is no
  `getPathnameList`.
- THE CATALOG RETURNS ONE PATHNAME PER STORAGE BLOCK, NOT ONE PER RECORD. A
  1HOUR series covering 01 Oct -> 01 May comes back eight times, once per month,
  each with its own D part. Reading a block-specific pathname returns ONLY that
  block. Blank the D part before reading (`parts[4] = ""`) and let DSS assemble
  the record; `read_ts` on `//B/C//1HOUR/F/` returns the whole thing. Indexing
  the catalog by (B, C, F) and keeping the first path seen returns roughly one
  month per record, and since catalog order differs per record, DIFFERENT
  months for different records -- two series from the same run land on
  non-overlapping dates and look like they used different mappings. This bit
  #Extract_Ensemble_To_Timeseries.py in Aug 2026; it now blanks the D part and
  stops if every member comes back under half its mapped hours.
- `import pydsstools` emits a `rasterio` ImportError on stderr. Harmless.
- Pathname case is inconsistent across records (`1Hour` vs `1HOUR`). Compare
  case-insensitively or you will report false mismatches.
- The header byte at offset 12 is `6` for a v6 file and `0` for a v7 file.
  Passing ANY explicit `version=` bypasses the auto-detect guard; the file is
  still opened with its own real version, so `version=7` on a v6 file works.
- DSS stamps are END of period. A 1Hour value at 01:00 covers hour 0 of that
  day; a 1Day value at midnight is the PREVIOUS calendar day. Get this wrong and
  everything is off by one.


## pandas 3 trap (silent wrong answers)

`DatetimeIndex.asi8` returns **microseconds** in pandas 3, while
`pd.Timestamp.value` returns **nanoseconds**. Mixing them in an `np.interp`
call is off by 1000x and does not raise -- it silently returns a flat line.
This produced a rule curve stuck at its minimum value until it was caught by
checking the output against the input anchor points.

Interpolate on unit-free floats instead:

    epoch = pd.Timestamp(1900, 1, 1)
    x = (index - epoch) / pd.Timedelta(hours=1)

The user's machine may be on pandas 2, where the mixed version happens to work.
Do not rely on that.

Always verify an interpolated series against its input anchors before shipping.


## Anything that runs on BOTH machines must be API-agnostic

Three pydsstools differences have now cost a round trip each: `ts.times`
(HecTime here, strings there), midnight-as-2400 date strings, and `search_path`
(present here, absent there). The pattern is the same every time -- code that
works in the sandbox fails on Windows because a pydsstools API differs.

Before shipping anything that touches pydsstools, ask which API is being called
and whether a different build might name it differently. Prefer `getattr` with
fallbacks over a direct call. There is no way to test the user's build from
here, so defensive code is the only protection.


## Repo conventions

- Shared modules live in `/Modules` at the repo root. Nothing is duplicated per
  project. Scripts add it to `sys.path` with a `".."` hop count matching their
  own depth.
- `REPO_ROOT` is derived from `__file__`, never hardcoded.
- Diagnostic output goes in `<project>/output/diagnostics/`, kept separate from
  `<project>/output/` which holds results that feed later steps.
- Scripts `os.chdir` to their own folder and use relative paths from there.
- `CAS_Unreg_FF/data/obsData.dss` is the canonical shared store for BOTH
  projects. Several scripts WRITE into it. Back it up before running those.
- User prefers flat functions, no classes, no argparse, hardcoded paths at the
  top, and full scripts rather than snippets.


## Ensemble workflow order (matters)

`#Create_ObsRC_Ensembles.py` times its windows off the REGULATED Castle Rock
peak, which only exists after the WCM_RC run has been extracted. The order is:

    1. #Create_Unreg_Ensembles.py        -> ensemble_unreg.dss   (WCM rule curve)
    2. run ResSim
    3. #Extract_Ensemble_To_Timeseries.py with SET_NAME="ResSim_WCM_RC"
    4. #Create_ObsRC_Ensembles.py        -> ensemble_obs_rc.dss  (observed pool)
    5. run ResSim
    6. #Extract_Ensemble_To_Timeseries.py with SET_NAME="ResSim_Obs_RC"

Step 4 fails loudly if step 3 has not been done. Do NOT flip
FALLBACK_TO_INFLOW_SUM to work around it -- that silently substitutes
UNREGULATED peak timing, which is a different analysis.


## What the two ensemble sets are for

Both runs use the SAME observed hydrology and the SAME current operating rules.
Only the reservoir's starting pool differs, which is what isolates the value of
the starting condition:

    ResSim_WCM_RC  pool starts at the WCM rule curve   (Oct -> May windows)
    ResSim_Obs_RC  pool starts at the OBSERVED elevation at event onset
                   (31-day windows on the event that made the WCM_RC peak)

The observed elevation in Obs_RC is a LOOKBACK record. ResSim reads only the
value at the simulation start from it and then lets the release rules control
the pool -- a full-length record does NOT pin the pool across the run. It is
written full length (ELEV_EXTENT="full") on purpose, so the simulation start
time can be shifted anywhere inside the window without rebuilding the ensemble.

The endpoint is three sets of Castle Rock regulated peaks to compare: the USGS
record, ResSim_WCM_RC, and ResSim_Obs_RC. The USGS annual peaks are already in
obsData.dss at
/COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW-ANNUAL PEAK//IR-CENTURY/USGS/.


## Obs_RC is limited to WY1974 onward

//MOS/ELEV//1DAY/USGS/ starts 02 Oct 1973, so 45 of the 98 water years have no
observed starting pool and cannot be part of the Obs_RC run. The script skips
them and lists them. That is a real limit of the record, not a bug -- Mossyrock
only came online in 1968 anyway.

Symptom if this check is removed: DSS accepts an all-missing record on write but
those records then fail to read back with "Error code -1". A record that is
entirely -901 is effectively corrupt. Never write one.


## The mismatched-run trap (already caught once)

Reassembling run A's simulation with run B's mapping does NOT error. It stamps
A's members onto B's dates and writes a series that looks entirely plausible --
right length, right shape, wrong dates. It was caught only because the
pass-through FLOW-IN record disagreed with ResSimInflows.dss.

SIM_DSS now lives in CONFIG_BY_SET alongside the mapping, keyed by SET_NAME, so
the pair cannot drift apart. A guard also stops any record whose members return
more than WINDOW_TOLERANCE x the mapped hours.

Diagnostic that found it: take a member's reassembled block, scan the source
record for the offset where it matches, and see where it really came from. If
that date is outside the member's own window, the wrong simulation was read.


## The synthetic ensemble is hand-edited after it is built

`ensemble_synthetic.dss` is NOT purely a script product. After
`#Create_Synthetic_Ensembles.py` writes it, the user chops a mini peak out of
the **Dec1977 and Nov1986** members by hand in DSSVue (both the
`MOSSYROCK/FLOW-IN` and `CASTLE ROCK/FLOW-LOCAL` records), then loads the
edited file into ResSim.

Why the bump is there: matching a peak AND a 5-day volume on a sharp event
forces `f_out` above `f_peak` (1.9x vs 1.1x for Dec1977), so a small rise the
observed hydrograph already had gets stretched into what looks like a second
flood. At the 500-year target it reached 0.43-0.45 of the peak.

**`#Create_Synthetic_Ensembles.py` deletes OUT_DSS before writing, so running
it destroys the chop.** Nothing downstream notices — the members still hit
their peak and volume targets with the bump present. Never re-run that script
without telling the user the hand edit has to be redone.

A scripted alternative was tried in Aug 2026 (raising the flow weight to a
power to concentrate the multiplier near the peak) and reverted: the amount of
damping that is safe — anything past ~1.4 lets a shoulder scale above the peak
and invents a new maximum — only moved the bump 1-2%, which the user judged
not worth the complexity. Do not re-propose it without a materially better
idea.

**The repo is in a mixed state and it is easy to get this wrong.** As of
commit `db5e090` the RESULTS are post-chop (`ResSim_Synth.dss`,
`synthetic_results.csv`) but the INPUT is not — `ensemble_synthetic.dss` is
still the unchopped script output, because the chopped ensemble was never
uploaded. Re-running ResSim from the repo copy silently reproduces the pre-chop
answer. Do not assume the two agree.

Check by hydrograph, never by commit date. Nov1986 500-yr (member 47, synthetic
WY1847) against its 228,861 cfs peak:

    ensemble_synthetic.dss   103,600 cfs  0.45   <- unchopped build
    ResSim_Synth.dss          97,746 cfs  0.43   <- post-chop (current)
                             105,586 cfs  0.46   <- what the pre-chop run read

The chop moved four members: Dec1977 250yr/500yr/beyond and Nov1986 beyond,
regulated peaks down 1.4% to 5.4%. The other 44 are identical. It works by
holding the pool lower going into the main peak, so the effect is on the
REGULATED side even where the unregulated hydrograph barely moved (Dec1977's
routed unreg bump is 0.38 of its peak either way).


## Open issue as of 7 Aug 2026

`ObsData_RegUnreg.dss` was deleted from `CAS_Reg_Unreg/data` as "duplicative",
but 644 of its 1,108 records never made it into `obsData.dss`:

    210 blocks  //MAY/FLOW-LOCAL-SHAPED//1Hour/CWMS/
    212 blocks  //MOS/FLOW-IN-CALC-CLEANED-VOLCOR//1Hour/CWMS/
    167 blocks  //MOS/ELEV-RULECURVE//1Hour/CENWP-CALC/
     54 blocks  //MOS/ELEV-USGS//1Day/USGS/
      1 block   //MOS/Location Info////

Only the Mayfield USGS gage carried over, and it was already present with wider
coverage. There are no `//MAY/` records of any kind in `obsData.dss`.

Nothing is lost -- the full v7 file is recoverable from git at `fd5909c~1`:

    git show fd5909c~1:CAS_Reg_Unreg/data/ObsData_RegUnreg.dss > old.dss

`CAS_Reg_Unreg/src/#MOS_Special_Release_MinFloodPool.py` still points at the
deleted file. It is deliberately NOT repointed: the volume-corrected inflows and
shaped Mayfield local it uses are being replaced by `ResSimInflows.dss` once the
ResSim inputs are final.
