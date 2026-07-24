# MOS_CASTLEROCK_PEAK_DATE_COMPARE.py

import os
import sys
import ssl
import certifi
import tempfile
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import dataretrieval.nwis as nwis
import cwms
from pydsstools.heclib.dss import HecDss


# =============================================================================
# USER SETTINGS
# =============================================================================
WRITE_CSV = False

START_DATE = "2008-01-01"
END_DATE = "2026-10-01"

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

UNREG_CSV_NAME = "../../Cowlitz_FF_DataPrep/data/2009_WY_MAX.csv"

DSS_FILE = r"../../Cowlitz_FF_DataPrep/data/results.dss"
UNREG_DSS_PATH = r"//CASTLEROCK_NWS/FLOW-UNREG/29Oct2008 - 28Mar2026/1Hour/UNREG_RESULTS/"

USGS_CASTLE_ROCK_SITE = "14243000"
USGS_FLOW_PARAMETER = "00060"

CWMS_TS_ID = "MOS.Stor.Inst.0.0.MIXED-COMPUTED-REV"

OUT_UNREG_CSV = "../diagnostics/Castlerock_unregulated_WY_max.csv"
OUT_DAILY_FLOW_CSV = "../diagnostics/CastleRock_USGS_daily_flow.csv"
OUT_IV_FLOW_CSV = "../diagnostics/CastleRock_USGS_IV_flow.csv"
OUT_USGS_PEAK_CSV = "../diagnostics/CastleRock_USGS_peak_flow_record.csv"
OUT_MOS_STORAGE_CSV = "../diagnostics/MOS_storage.csv"
OUT_PEAK_DATE_COMPARE_CSV = (
    "../diagnostics/Castlerock_unreg_dss_vs_usgs_regulated_peak_date_compare_WY1988_present.csv"
)
# =============================================================================


if getattr(sys, "frozen", False):
    BASE_DIR = Path(sys.argv[0]).resolve().parent
else:
    BASE_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------------
# SSL certificate setup
# ---------------------------------------------------------------------
pem_path = os.path.join(tempfile.gettempdir(), "corp_plus_certifi.pem")


def build_windows_ca_bundle(target_pem):
    base_bundle = certifi.where()

    with open(base_bundle, "rb") as src, open(target_pem, "wb") as dst:
        dst.write(src.read())

        try:
            for cert_tuple in ssl.enum_certificates("ROOT"):
                der_bytes = cert_tuple[0]
                pem_str = ssl.DER_cert_to_PEM_cert(der_bytes)
                dst.write(pem_str.encode("ascii"))
        except Exception as e:
            print(f"[WARNING] Error reading Windows ROOT store: {e}")

    return target_pem


if not os.path.exists(pem_path):
    try:
        bundle_path = build_windows_ca_bundle(pem_path)
    except Exception:
        bundle_path = certifi.where()
else:
    bundle_path = pem_path

os.environ["REQUESTS_CA_BUNDLE"] = bundle_path
print(f"[INFO] Using CA bundle: {bundle_path}")


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def water_year(dt):
    dt = pd.to_datetime(dt)
    return dt.year + 1 if dt.month >= 10 else dt.year


def make_tz_naive_index(obj):
    obj = obj.copy()
    obj.index = pd.to_datetime(obj.index)

    if getattr(obj.index, "tz", None) is not None:
        obj.index = obj.index.tz_convert("UTC").tz_localize(None)

    return obj.sort_index()


def clean_numeric_columns(df, columns):
    df = df.copy()

    for col in columns:
        df[col] = (
            df[col]
            .astype(str)
            .str.replace(",", "", regex=False)
            .str.strip()
        )
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def extract_usgs_flow(df, parameter_cd, flow_name):
    if df.empty:
        return pd.Series(dtype=float, name=flow_name)

    df = make_tz_naive_index(df)

    flow_cols = [
        c for c in df.columns
        if str(c).startswith(parameter_cd)
        and not str(c).endswith("_cd")
    ]

    if not flow_cols:
        raise ValueError(f"No USGS flow column found. Columns: {list(df.columns)}")

    s = pd.to_numeric(df[flow_cols[0]], errors="coerce")
    s = s.replace([-999999, -99999, -9999, -902, -901], np.nan)
    s = s.dropna()
    s = s[s > 0]
    s.name = flow_name

    return s


# ---------------------------------------------------------------------
# Import unregulated WY max CSV
# ---------------------------------------------------------------------
unreg_path = BASE_DIR / UNREG_CSV_NAME

unreg_max = pd.read_csv(unreg_path)
unreg_max.columns = [c.strip() for c in unreg_max.columns]

expected_cols = ["WY", "Peak", "One_day", "Three_Day", "Five_Day"]
missing_cols = [c for c in expected_cols if c not in unreg_max.columns]

if missing_cols:
    raise ValueError(f"Missing expected columns in {UNREG_CSV_NAME}: {missing_cols}")

unreg_max = clean_numeric_columns(
    unreg_max,
    ["WY", "Peak", "One_day", "Three_Day", "Five_Day"],
)

unreg_max["WY"] = unreg_max["WY"].astype("Int64")

print("\nImported unregulated WY max CSV:")
print(f"  File: {unreg_path}")
print(f"  Records: {len(unreg_max)}")
print(f"  WY range: {unreg_max['WY'].min()} to {unreg_max['WY'].max()}")


# ---------------------------------------------------------------------
# Import unregulated hourly Castlerock DSS record
# ---------------------------------------------------------------------
def import_dss_timeseries(dss_file, pathname, value_name):
    print("\nImporting DSS time series:")
    print(f"  File: {dss_file}")
    print(f"  Path: {pathname}")

    dss = HecDss.Open(dss_file)

    try:
        ts = dss.read_ts(pathname)

        df = pd.DataFrame({
            "DateTime": pd.to_datetime(ts.pytimes),
            value_name: pd.to_numeric(ts.values, errors="coerce"),
        })
    finally:
        dss.close()

    df = df.dropna()
    df = df[df[value_name] > -9000]
    df = df.set_index("DateTime").sort_index()

    print(f"  Records: {len(df)}")
    print(f"  First: {df.index.min()}")
    print(f"  Last:  {df.index.max()}")

    return df


unreg_dss = import_dss_timeseries(
    DSS_FILE,
    UNREG_DSS_PATH,
    "Unregulated_Flow_cfs",
)

unreg_dss["WY"] = unreg_dss.index.map(water_year)

unreg_peak_by_wy = (
    unreg_dss.loc[unreg_dss.groupby("WY")["Unregulated_Flow_cfs"].idxmax()]
    .reset_index()
    .rename(columns={
        "DateTime": "Unregulated_Peak_DateTime",
        "Unregulated_Flow_cfs": "Unregulated_Peak_cfs",
    })
)


# ---------------------------------------------------------------------
# Download USGS daily Castle Rock flow
# ---------------------------------------------------------------------
print("\nDownloading USGS daily Castle Rock flow:")
print(f"  Site: {USGS_CASTLE_ROCK_SITE}")
print(f"  {START_DATE} to {END_DATE}")

daily_raw = nwis.get_record(
    sites=USGS_CASTLE_ROCK_SITE,
    service="dv",
    start=START_DATE,
    end=END_DATE,
    parameterCd=USGS_FLOW_PARAMETER,
)

daily_flow = extract_usgs_flow(
    daily_raw,
    USGS_FLOW_PARAMETER,
    "CastleRock_Daily_Flow_cfs",
)

print(f"  Downloaded {len(daily_flow)} daily records.")
if not daily_flow.empty:
    print(f"  First: {daily_flow.index.min()}")
    print(f"  Last:  {daily_flow.index.max()}")


# ---------------------------------------------------------------------
# Download USGS instantaneous Castle Rock flow
# ---------------------------------------------------------------------
print("\nDownloading USGS instantaneous Castle Rock flow:")
print(f"  Site: {USGS_CASTLE_ROCK_SITE}")
print(f"  {START_DATE} to {END_DATE}")

iv_raw = nwis.get_record(
    sites=USGS_CASTLE_ROCK_SITE,
    service="iv",
    start=START_DATE,
    end=END_DATE,
    parameterCd=USGS_FLOW_PARAMETER,
)

iv_flow = extract_usgs_flow(
    iv_raw,
    USGS_FLOW_PARAMETER,
    "CastleRock_IV_Flow_cfs",
)

print(f"  Downloaded {len(iv_flow)} IV records.")
if not iv_flow.empty:
    print(f"  First: {iv_flow.index.min()}")
    print(f"  Last:  {iv_flow.index.max()}")


# ---------------------------------------------------------------------
# Download USGS regulated annual peak flow record
# ---------------------------------------------------------------------
def download_usgs_peak_record(site):
    print("\nDownloading USGS regulated annual peak flow record:")
    print(f"  Site: {site}")

    url = "https://nwis.waterdata.usgs.gov/nwis/peak"

    params = {
        "site_no": site,
        "agency_cd": "USGS",
        "format": "rdb",
    }

    response = requests.get(
        url,
        params=params,
        verify=os.environ.get("REQUESTS_CA_BUNDLE", True),
        timeout=60,
    )
    response.raise_for_status()

    df = pd.read_csv(
        StringIO(response.text),
        sep="\t",
        comment="#",
        dtype=str,
    )

    df = df[df["peak_dt"] != "10d"].copy()

    df["USGS_Regulated_Peak_Date"] = pd.to_datetime(
        df["peak_dt"],
        errors="coerce",
    )

    df["USGS_Regulated_Peak_cfs"] = pd.to_numeric(
        df["peak_va"],
        errors="coerce",
    )

    df = df.dropna(
        subset=[
            "USGS_Regulated_Peak_Date",
            "USGS_Regulated_Peak_cfs",
        ]
    )

    df["WY"] = df["USGS_Regulated_Peak_Date"].map(water_year)

    print(f"  Downloaded {len(df)} annual peaks.")
    print(f"  First WY: {df['WY'].min()}")
    print(f"  Last WY:  {df['WY'].max()}")

    return df


peak_record = download_usgs_peak_record(USGS_CASTLE_ROCK_SITE)

usgs_peak_by_wy = peak_record[
    [
        "WY",
        "USGS_Regulated_Peak_Date",
        "USGS_Regulated_Peak_cfs",
    ]
].copy()


# ---------------------------------------------------------------------
# Download CWMS Mossyrock storage
# ---------------------------------------------------------------------
def download_cwms_timeseries(ts_id, start_date, end_date, office="nws"):
    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)

    api_root = f"https://wm.{office}.ds.usace.army.mil:8243/nwdp-data/"
    cwms.api.init_session(api_root=api_root)

    print("\nDownloading CWMS Mossyrock storage:")
    print(f"  TS ID: {ts_id}")
    print(f"  {start_dt} to {end_dt}")

    data = cwms.get_timeseries(
        ts_id,
        office_id="NWDP",
        begin=start_dt,
        end=end_dt,
    ).df

    if data.empty:
        print("[WARNING] CWMS storage data is empty.")
        return pd.DataFrame(columns=["MOS_Storage_acre_ft"])

    if "date-time" in data.columns:
        data = data.set_index("date-time")

    data = data.iloc[:, [0]].copy()
    data.columns = ["MOS_Storage_acre_ft"]

    data["MOS_Storage_acre_ft"] = pd.to_numeric(
        data["MOS_Storage_acre_ft"],
        errors="coerce",
    )

    data = data.dropna()
    data = data[data["MOS_Storage_acre_ft"] > -9000]
    data = make_tz_naive_index(data)

    print(f"  Downloaded {len(data)} storage records.")
    print(f"  First: {data.index.min()}")
    print(f"  Last:  {data.index.max()}")

    return data


mos_storage = download_cwms_timeseries(
    CWMS_TS_ID,
    START_DATE,
    END_DATE,
)


# ---------------------------------------------------------------------
# Compare unregulated DSS peak dates to USGS regulated peak dates
# Calculate flow reduction and Mossyrock storage changes
# Regress flow reduction against 1-, 3-, and 5-day storage increases
# ---------------------------------------------------------------------
peak_date_compare = unreg_peak_by_wy.merge(
    usgs_peak_by_wy,
    on="WY",
    how="left",
)

peak_date_compare = peak_date_compare[
    peak_date_compare["WY"] >= 2008
].copy()

peak_date_compare["Unregulated_Peak_Date"] = (
    peak_date_compare["Unregulated_Peak_DateTime"].dt.floor("D")
)

peak_date_compare["Date_Difference_Days"] = (
    peak_date_compare["Unregulated_Peak_Date"]
    - peak_date_compare["USGS_Regulated_Peak_Date"]
).dt.days

peak_date_compare["Peak_Difference_cfs"] = (
    peak_date_compare["Unregulated_Peak_cfs"]
    - peak_date_compare["USGS_Regulated_Peak_cfs"]
)

peak_date_compare["Peak_Percent_Difference"] = (
    peak_date_compare["Peak_Difference_cfs"]
    / peak_date_compare["USGS_Regulated_Peak_cfs"]
    * 100.0
)


# ---------------------------------------------------------------------
# Regulated daily flow at Castle Rock on each unregulated peak date
# ---------------------------------------------------------------------
regulated_daily = daily_flow.copy()
regulated_daily.index = pd.to_datetime(regulated_daily.index).floor("D")
regulated_daily = regulated_daily.groupby(regulated_daily.index).mean()
regulated_daily.name = "Regulated_Flow_On_Unreg_Peak_Date_cfs"

peak_date_compare = peak_date_compare.merge(
    regulated_daily,
    left_on="Unregulated_Peak_Date",
    right_index=True,
    how="left",
)

peak_date_compare["Flow_Reduction_cfs"] = (
    peak_date_compare["Unregulated_Peak_cfs"]
    - peak_date_compare["Regulated_Flow_On_Unreg_Peak_Date_cfs"]
)


# ---------------------------------------------------------------------
# Mossyrock storage increase over 1, 3, and 5 days before unregulated peak
# ---------------------------------------------------------------------
storage = mos_storage.copy()
storage = storage.sort_index()

storage_col = "MOS_Storage_acre_ft"

def storage_asof(dt):
    dt = pd.to_datetime(dt)
    prior = storage.loc[storage.index <= dt, storage_col]

    if prior.empty:
        return np.nan

    return prior.iloc[-1]


for n_days in [1, 3, 5]:
    peak_date_compare[f"MOS_Storage_at_Peak_acft"] = (
        peak_date_compare["Unregulated_Peak_DateTime"].apply(storage_asof)
    )

    peak_date_compare[f"MOS_Storage_{n_days}Day_Before_acft"] = (
        peak_date_compare["Unregulated_Peak_DateTime"]
        .apply(lambda x: storage_asof(pd.to_datetime(x) - pd.Timedelta(days=n_days)))
    )

    peak_date_compare[f"MOS_Storage_Increase_{n_days}Day_acft"] = (
        peak_date_compare["MOS_Storage_at_Peak_acft"]
        - peak_date_compare[f"MOS_Storage_{n_days}Day_Before_acft"]
    )


# ---------------------------------------------------------------------
# Simple OLS regression helper
# y = intercept + slope * x
# ---------------------------------------------------------------------
def simple_ols(df, x_col, y_col):
    reg_df = df[[x_col, y_col]].dropna().copy()

    if len(reg_df) < 3:
        return {
            "Predictor": x_col,
            "N": len(reg_df),
            "Intercept": np.nan,
            "Slope": np.nan,
            "R2": np.nan,
            "RMSE_cfs": np.nan,
        }

    x = reg_df[x_col].to_numpy(dtype=float)
    y = reg_df[y_col].to_numpy(dtype=float)

    X = np.column_stack([np.ones(len(x)), x])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]

    y_hat = X @ beta
    residuals = y - y_hat

    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = np.sqrt(np.mean(residuals ** 2))

    return {
        "Predictor": x_col,
        "N": len(reg_df),
        "Intercept": beta[0],
        "Slope": beta[1],
        "R2": r2,
        "RMSE_cfs": rmse,
    }


regression_results = pd.DataFrame([
    simple_ols(
        peak_date_compare,
        "MOS_Storage_Increase_1Day_acft",
        "Flow_Reduction_cfs",
    ),
    simple_ols(
        peak_date_compare,
        "MOS_Storage_Increase_3Day_acft",
        "Flow_Reduction_cfs",
    ),
    simple_ols(
        peak_date_compare,
        "MOS_Storage_Increase_5Day_acft",
        "Flow_Reduction_cfs",
    ),
])

best_predictor = regression_results.sort_values(
    ["R2", "RMSE_cfs"],
    ascending=[False, True],
).iloc[0]


# ---------------------------------------------------------------------
# Optional combined regression:
# y = intercept + b1*1day + b3*3day + b5*5day
# ---------------------------------------------------------------------
combined_cols = [
    "MOS_Storage_Increase_1Day_acft",
    "MOS_Storage_Increase_3Day_acft",
    "MOS_Storage_Increase_5Day_acft",
    "Flow_Reduction_cfs",
]

combined_df = peak_date_compare[combined_cols].dropna().copy()

if len(combined_df) >= 5:
    y = combined_df["Flow_Reduction_cfs"].to_numpy(dtype=float)

    X = combined_df[
        [
            "MOS_Storage_Increase_1Day_acft",
            "MOS_Storage_Increase_3Day_acft",
            "MOS_Storage_Increase_5Day_acft",
        ]
    ].to_numpy(dtype=float)

    X = np.column_stack([np.ones(len(X)), X])

    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_hat = X @ beta
    residuals = y - y_hat

    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    combined_regression = pd.DataFrame([{
        "N": len(combined_df),
        "Intercept": beta[0],
        "Slope_1Day": beta[1],
        "Slope_3Day": beta[2],
        "Slope_5Day": beta[3],
        "R2": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        "RMSE_cfs": np.sqrt(np.mean(residuals ** 2)),
    }])
else:
    combined_regression = pd.DataFrame()


peak_date_compare = peak_date_compare[
    [
        "WY",
        "Unregulated_Peak_DateTime",
        "Unregulated_Peak_Date",
        "Unregulated_Peak_cfs",
        "USGS_Regulated_Peak_Date",
        "USGS_Regulated_Peak_cfs",
        "Date_Difference_Days",
        "Regulated_Flow_On_Unreg_Peak_Date_cfs",
        "Flow_Reduction_cfs",
        "MOS_Storage_at_Peak_acft",
        "MOS_Storage_1Day_Before_acft",
        "MOS_Storage_Increase_1Day_acft",
        "MOS_Storage_3Day_Before_acft",
        "MOS_Storage_Increase_3Day_acft",
        "MOS_Storage_5Day_Before_acft",
        "MOS_Storage_Increase_5Day_acft",
        "Peak_Difference_cfs",
        "Peak_Percent_Difference",
    ]
].sort_values("WY")


print("\nUnregulated peak, regulated flow, and Mossyrock storage changes:")
print(peak_date_compare.to_string(index=False))

print("\nSimple regression results:")
print(regression_results.to_string(index=False))

print("\nBest single predictor:")
print(best_predictor.to_string())

print("\nCombined regression results:")
if combined_regression.empty:
    print("Not enough complete records for combined regression.")
else:
    print(combined_regression.to_string(index=False))

# ---------------------------------------------------------------------
# Optional CSV writes
# ---------------------------------------------------------------------
if WRITE_CSV:
    unreg_max.to_csv(BASE_DIR / OUT_UNREG_CSV, index=False)
    daily_flow.to_csv(BASE_DIR / OUT_DAILY_FLOW_CSV)
    iv_flow.to_csv(BASE_DIR / OUT_IV_FLOW_CSV)
    peak_record.to_csv(BASE_DIR / OUT_USGS_PEAK_CSV, index=False)
    mos_storage.to_csv(BASE_DIR / OUT_MOS_STORAGE_CSV)
    peak_date_compare.to_csv(BASE_DIR / OUT_PEAK_DATE_COMPARE_CSV, index=False)

    print("\nWrote:")
    print(f"  {BASE_DIR / OUT_UNREG_CSV}")
    print(f"  {BASE_DIR / OUT_DAILY_FLOW_CSV}")
    print(f"  {BASE_DIR / OUT_IV_FLOW_CSV}")
    print(f"  {BASE_DIR / OUT_USGS_PEAK_CSV}")
    print(f"  {BASE_DIR / OUT_MOS_STORAGE_CSV}")
    print(f"  {BASE_DIR / OUT_PEAK_DATE_COMPARE_CSV}")
else:
    print("\nCSV writing is OFF. Set WRITE_CSV = True to write output files.")