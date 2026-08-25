# ecology_io.py
# -*- coding: utf-8 -*-
"""Reading WA Dept of Ecology gage exports, with the quality column honoured.

WHY THIS IS A SHARED MODULE
    The row regex used to live as a copy-pasted constant in four Coweeman
    scripts, and all four copies carried the same bug for months (see below).
    That is exactly the failure mode /Modules exists to prevent, so the parser
    lives here now and nothing re-declares it.

THE BUG THIS MODULE EXISTS TO STOP RECURRING
    Ecology's exports are FIXED WIDTH, and the value column is left BLANK
    whenever no value is reported -- the quality code then sits alone on the
    line:

        10/01/2016   00:00              35.2         2    value 35.2, code 2
        03/16/2017   05:00                          254    NO VALUE, code 254

    The old pattern made the value mandatory and the quality optional:

        r"...(-?[\\d.]+)(?:\\s+(\\S+))?\\s*$"          # WRONG

    so the second form parsed as a discharge of 254.0 cfs. Code 254 is
    "Rating table exceeded (data will not be reported)", so on gage 26C075 --
    whose rating tops out near 3,400 cfs -- the HIGHEST flows in the record
    silently became one of its LOWEST values. Nothing raised, nothing logged.

    The pattern below inverts it: the value is optional, the quality code is a
    REQUIRED integer. Given "   254" the engine first tries value="254", finds
    no quality token after it, and backtracks to value=None, quality=254 --
    the correct reading. Verified against all 15 cached 26C075 files: 462,445
    data rows parse and not one lacks a quality code.

CODES THAT CARRY NO USABLE VALUE  (MISSING_CODES -> NaN, never a number)
    254  Rating table exceeded (data will not be reported). Value column blank.
    151  Data Missing. These rows DO carry a number and it decays to 0.0 as
         the sensor fails -- on 09 Dec 2015 that puts a 0 cfs Coweeman in the
         middle of the second largest Cowlitz event in the record.

    Both are NaN here. A censored reading is not a low flow; on this gage it
    is an UNKNOWN AND HIGH flow, and treating it as a number biases every
    coincident ratio toward zero on precisely the largest events.

CODES THAT ARE A REAL ATTEMPT AT A VALUE BUT NOT A GAGED MEASUREMENT
    10 above rating (within 2x), 50 estimated, 77 estimated from another
    station, 82 interpolated across a gap, 100 modeled, 140 not yet checked,
    and 160 -- which appears in NO legend in ANY cached file yet carries 1,012
    readings spanning 2,450-8,020 cfs, including the record maximum. Callers
    that care should test against TRUSTED_CODES rather than assume.
"""

import glob
import os
import re

import numpy as np
import pandas as pd

# Value optional, quality REQUIRED -- see module docstring.
ECOLOGY_ROW_DT = re.compile(
    r"^\s*(\d{1,2}/\d{1,2}/\d{4})\s+(\d{1,2}:\d{2})\s+(-?[\d.]+)?\s+(\d+)\s*$")
ECOLOGY_ROW_D = re.compile(
    r"^\s*(\d{1,2}/\d{1,2}/\d{4})\s+(-?[\d.]+)?\s+(\d+)\s*$")

MISSING_CODES = (151, 254)
TRUSTED_CODES = (1, 2, 3, 8)

CODE_MEANING = {
    1: "good, reviewed",
    2: "good, provisional",
    3: "good, provisional - edited",
    8: "below rating",
    10: "ABOVE RATING (within 2x)",
    50: "estimated",
    77: "estimated from another station",
    82: "interpolated across gap",
    100: "MODELED",
    140: "not yet checked",
    151: "DATA MISSING",
    160: "UNDOCUMENTED CODE (in no legend in any file)",
    254: "RATING EXCEEDED, not reported",
}


def parse_ecology_text(text):
    """(value, quality) Series from one Ecology export.

    Handles the 15-minute form (FM, DATE + TIME) and the mean-daily form
    (DV, DATE only). Values under MISSING_CODES come back NaN with the code
    preserved, so a caller can tell "unknown and above rating" from "absent".
    """
    stamps, values, quality = [], [], []
    for line in text.splitlines():
        line = line.rstrip("\r")
        match = ECOLOGY_ROW_DT.match(line)
        if match:
            date, clock, value, code = match.groups()
            stamps.append("%s %s" % (date, clock))
        else:
            match = ECOLOGY_ROW_D.match(line)
            if not match:
                continue
            date, value, code = match.groups()
            stamps.append(date)
        values.append(np.nan if value is None else float(value))
        quality.append(int(code))
    if len(values) < 2:
        return None, None

    # Kept as one frame so the code stays with its value POSITIONALLY.
    # Aligning the two afterwards by timestamp breaks the moment a file
    # repeats one, which a daily file does trivially.
    frame = pd.DataFrame({"stamp": stamps, "value": values, "qual": quality})
    frame["stamp"] = pd.to_datetime(frame["stamp"], format="mixed",
                                    errors="coerce")
    frame = frame.dropna(subset=["stamp"])
    frame = frame.drop_duplicates(subset="stamp", keep="first")
    frame = frame.sort_values("stamp")

    # Ecology also writes missing record as a large negative rather than a
    # blank, on top of the quality codes.
    frame.loc[frame["value"] <= -900.0, "value"] = np.nan
    frame.loc[frame["qual"].isin(MISSING_CODES), "value"] = np.nan

    index = pd.DatetimeIndex(frame["stamp"])
    return (pd.Series(frame["value"].to_numpy(), index=index),
            pd.Series(frame["qual"].to_numpy(), index=index))


def read_ecology_cache(cache_dir, pattern="*_FM.txt", label="Ecology gage",
                       verbose=True):
    """(value, quality) for every cached export matching pattern."""
    files = sorted(glob.glob(os.path.join(cache_dir, pattern)))
    if not files:
        raise SystemExit(
            "No Ecology files matching %s in %s.\nRun #Coweeman_Timing.py "
            "first; it downloads and caches them."
            % (pattern, os.path.abspath(cache_dir)))
    values, quals = [], []
    for path in files:
        with open(path, "r", errors="replace") as handle:
            value, qual = parse_ecology_text(handle.read())
        if value is not None:
            values.append(value)
            quals.append(qual)
    value = pd.concat(values)
    qual = pd.concat(quals)
    keep = ~value.index.duplicated(keep="first")
    value, qual = value[keep].sort_index(), qual[keep].sort_index()

    if verbose:
        missing = qual.isin(MISSING_CODES)
        qualified = (~qual.isin(TRUSTED_CODES)) & (~missing)
        print("   %-10s: %d readings from %d file(s), %s to %s"
              % (label, len(value), len(files), value.index.min().date(),
                 value.index.max().date()))
        print("               %d over %d day(s) carry NO usable value (codes "
              "%s) and are held as unknown, not as a number"
              % (int(missing.sum()), qual.index[missing].normalize().nunique(),
                 ", ".join(str(c) for c in MISSING_CODES)))
        print("               %d are estimated/modeled/above-rating rather "
              "than gaged in rating" % int(qualified.sum()))
    return value, qual


def censored_spans(qual, lo=None, hi=None):
    """Contiguous runs carrying no usable value, as (start, end) stamps."""
    flag = qual.isin(MISSING_CODES)
    if lo is not None:
        flag = flag.loc[lo:hi]
    if not flag.any():
        return []
    stamps = flag.index[flag]
    breaks = np.flatnonzero(np.diff(stamps.values) > np.timedelta64(20, "m"))
    starts = np.r_[0, breaks + 1]
    ends = np.r_[breaks, len(stamps) - 1]
    return [(stamps[a], stamps[b]) for a, b in zip(starts, ends)]


def resample_censor_aware(value, qual, rule="1h", how="mean"):
    """Resample, but return NaN for any bin holding a censored reading.

    A plain .resample().mean() skips NaN, so an hour where the flow climbed
    ABOVE the rating part-way through would come back as the average of only
    the in-rating readings that preceded it -- a number materially lower than
    the truth, produced at exactly the moment the tributary was largest. That
    is the same failure the parser bug caused, one level up, so it is blocked
    here too: if any reading in the bin was censored, the bin is unknown.
    """
    out = getattr(value.resample(rule), how)()
    blocked = qual.isin(MISSING_CODES).resample(rule).max().astype(bool)
    out[blocked.reindex(out.index, fill_value=False)] = np.nan
    return out
