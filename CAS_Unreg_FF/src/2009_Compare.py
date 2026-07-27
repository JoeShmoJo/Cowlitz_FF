from pathlib import Path

import pandas as pd
import plotly.graph_objects as go


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

PAST_RESULTS_CSV = PROJECT_DIR / "data" / "CAS_2009_WY_Max.csv"
CURRENT_RESULTS_CSV = PROJECT_DIR / "output" / "wy_record_ssp.csv"
OUTPUT_HTML = PROJECT_DIR / "output" / "CAS_2009_comparison.html"


DURATIONS = {
    "Peak": {
        "past_column": "Peak",
        "current_column": "Peak",
        "label": "Peak",
    },
    "One_day": {
        "past_column": "1",
        "current_column": "One_day",
        "label": "1-Day",
    },
    "Three_Day": {
        "past_column": "3",
        "current_column": "Three_Day",
        "label": "3-Day",
    },
    "Five_Day": {
        "past_column": "5",
        "current_column": "Five_Day",
        "label": "5-Day",
    },
}


def read_csv_flexible(csv_path: Path) -> pd.DataFrame:
    """
    Read a comma-delimited or tab-delimited CSV.

    The delimiter is determined from the header so commas within numeric
    values such as "41,100" do not cause incorrect delimiter detection.
    """
    if not csv_path.exists():
        raise FileNotFoundError(f"Input file not found: {csv_path}")

    with csv_path.open("r", encoding="utf-8-sig") as file:
        header = file.readline()

    delimiter = "\t" if "\t" in header else ","

    return pd.read_csv(
        csv_path,
        sep=delimiter,
        dtype=str,
        encoding="utf-8-sig",
    )


def clean_numeric(series: pd.Series) -> pd.Series:
    """
    Convert numeric strings to numbers.

    Examples:
        "41,100" -> 41100
        ""       -> NaN
    """
    return pd.to_numeric(
        series.astype("string")
        .str.strip()
        .str.replace(",", "", regex=False),
        errors="coerce",
    )


def load_past_results() -> pd.DataFrame:
    df = read_csv_flexible(PAST_RESULTS_CSV)
    df.columns = df.columns.str.strip()

    required_columns = ["WY", "Peak", "1", "3", "5"]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"{PAST_RESULTS_CSV.name} is missing columns: "
            f"{missing_columns}"
        )

    df["WY"] = clean_numeric(df["WY"])
    df = df.dropna(subset=["WY"]).copy()
    df["WY"] = df["WY"].astype(int)

    for column in ["Peak", "1", "3", "5"]:
        df[column] = clean_numeric(df[column])

    return df.sort_values("WY").reset_index(drop=True)


def load_current_results() -> pd.DataFrame:
    df = read_csv_flexible(CURRENT_RESULTS_CSV)
    df.columns = df.columns.str.strip()

    required_columns = [
        "WY",
        "Peak",
        "One_day",
        "Three_Day",
        "Five_Day",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"{CURRENT_RESULTS_CSV.name} is missing columns: "
            f"{missing_columns}"
        )

    df["WY"] = clean_numeric(df["WY"])
    df = df.dropna(subset=["WY"]).copy()
    df["WY"] = df["WY"].astype(int)

    for column in [
        "Peak",
        "One_day",
        "Three_Day",
        "Five_Day",
    ]:
        df[column] = clean_numeric(df[column])

    return df.sort_values("WY").reset_index(drop=True)


def build_comparison_plot(
    past_df: pd.DataFrame,
    current_df: pd.DataFrame,
) -> go.Figure:
    figure = go.Figure()

    trace_duration_keys = []

    for duration_key, duration_info in DURATIONS.items():
        past_column = duration_info["past_column"]
        current_column = duration_info["current_column"]
        duration_label = duration_info["label"]

        past_data = (
            past_df[["WY", past_column]]
            .dropna(subset=[past_column])
            .copy()
        )

        current_data = (
            current_df[["WY", current_column]]
            .dropna(subset=[current_column])
            .copy()
        )

        figure.add_trace(
            go.Scatter(
                x=past_data["WY"].to_numpy(dtype=int),
                y=past_data[past_column].to_numpy(dtype=float),
                mode="markers",
                name=f"2009 Study — {duration_label}",
                legendgroup=duration_key,
                marker=dict(
                    symbol="circle",
                    size=6,
                    color="blue",
                    line=dict(width=0),
                ),
                hovertemplate=(
                    "<b>2009 Study</b><br>"
                    f"Duration: {duration_label}<br>"
                    "WY: %{x:.0f}<br>"
                    "Flow: %{y:,.0f} cfs"
                    "<extra></extra>"
                ),
            )
        )

        trace_duration_keys.append(duration_key)

        figure.add_trace(
            go.Scatter(
                x=current_data["WY"].to_numpy(dtype=int),
                y=current_data[current_column].to_numpy(dtype=float),
                mode="markers",
                name=f"Current Study — {duration_label}",
                legendgroup=duration_key,
                marker=dict(
                    symbol="circle",
                    size=10,
                    color="rgba(0,0,0,0)",
                    line=dict(
                        color="red",
                        width=2,
                    ),
                ),
                hovertemplate=(
                    "<b>Current Study</b><br>"
                    f"Duration: {duration_label}<br>"
                    "WY: %{x:.0f}<br>"
                    "Flow: %{y:,.0f} cfs"
                    "<extra></extra>"
                ),
            )
        )

        trace_duration_keys.append(duration_key)

    buttons = [
        {
            "label": "All Durations",
            "method": "update",
            "args": [
                {
                    "visible": [True] * len(figure.data),
                },
                {
                    "title": (
                        "Current Study vs. 2009 Study — "
                        "All Durations"
                    ),
                },
            ],
        }
    ]

    for duration_key, duration_info in DURATIONS.items():
        visible = [
            trace_key == duration_key
            for trace_key in trace_duration_keys
        ]

        buttons.append(
            {
                "label": duration_info["label"],
                "method": "update",
                "args": [
                    {
                        "visible": visible,
                    },
                    {
                        "title": (
                            "Current Study vs. 2009 Study — "
                            f"{duration_info['label']}"
                        ),
                    },
                ],
            }
        )

    all_water_years = pd.concat(
        [
            past_df["WY"],
            current_df["WY"],
        ],
        ignore_index=True,
    )

    min_water_year = int(all_water_years.min())
    max_water_year = int(all_water_years.max())

    figure.update_layout(
        title="Current Study vs. 2009 Study — All Durations",
        xaxis_title="Water Year",
        yaxis_title="Flow, cfs",
        hovermode="x unified",
        template="plotly_white",
        legend={
            "title": "Study and Duration",
            "groupclick": "toggleitem",
        },
        updatemenus=[
            {
                "type": "buttons",
                "direction": "right",
                "buttons": buttons,
                "x": 0.0,
                "xanchor": "left",
                "y": 1.14,
                "yanchor": "top",
                "showactive": True,
            }
        ],
        margin={
            "t": 130,
            "r": 40,
            "b": 70,
            "l": 90,
        },
    )

    figure.update_xaxes(
        type="linear",
        range=[
            min_water_year - 1,
            max_water_year + 1,
        ],
        tickmode="linear",
        tick0=(min_water_year // 5) * 5,
        dtick=5,
        tickformat="d",
        showgrid=True,
    )

    figure.update_yaxes(
        tickformat=",",
        rangemode="tozero",
        showgrid=True,
    )

    return figure


def main() -> None:
    past_df = load_past_results()
    current_df = load_current_results()

    print(
        "Past study water years:",
        past_df["WY"].min(),
        "to",
        past_df["WY"].max(),
    )

    print(
        "Current study water years:",
        current_df["WY"].min(),
        "to",
        current_df["WY"].max(),
    )

    figure = build_comparison_plot(
        past_df=past_df,
        current_df=current_df,
    )

    OUTPUT_HTML.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.write_html(
        OUTPUT_HTML,
        include_plotlyjs="cdn",
        full_html=True,
        auto_open=True,
    )

    print(f"Past-study records: {len(past_df):,}")
    print(f"Current-study records: {len(current_df):,}")
    print(f"Plot written to: {OUTPUT_HTML}")


if __name__ == "__main__":
    main()