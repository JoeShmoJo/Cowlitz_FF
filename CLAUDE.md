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
- `ts.times` yields `HecTime` objects. `.datetime()` is a **method** -- call it.
  Passing the objects to `pd.to_datetime` fails.
- Missing values: mask with `np.array(ts.nodata, dtype=bool)`, and also treat
  `<= -900` as missing. Write missing as `-3.4028234663852886e38`.
- `search_path('/*/*/*/*/*/*/')` is the catalog call. There is no
  `getPathnameList`.
- `import pydsstools` emits a `rasterio` ImportError on stderr. Harmless.
- Pathname case is inconsistent across records (`1Hour` vs `1HOUR`). Compare
  case-insensitively or you will report false mismatches.
- DSS stamps are END of period. A 1Hour value at 01:00 covers hour 0 of that
  day; a 1Day value at midnight is the PREVIOUS calendar day. Get this wrong and
  everything is off by one.


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
