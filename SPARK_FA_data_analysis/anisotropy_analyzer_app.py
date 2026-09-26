"""Fluorescence-anisotropy time-series analyzer (stand-alone desktop application).

Reads an Excel workbook of plate-reader time-series data

    column 1       : time (hours)
    column 2       : temperature (deg C)
    column 3 ...   : one column per well (e.g. fluorescence anisotropy, mA);
                     the first row holds the well names (D4, D5, ...)

and performs plotting, detrending, peak/trough detection, actograms,
interactive review of the peak/trough times (manual removal/addition),
period/phase regression on the peak/trough times and (damped) sine-curve
fitting for period, amplitude and phase.

The number of rows (duration) and of wells (columns) may differ between files.

Run with:   python anisotropy_analyzer_app.py
Requires:   numpy, pandas, scipy, matplotlib, openpyxl  (Tkinter ships with Python)

Layout of this file:
  1. Processing core  (no GUI code; can be used from scripts)
  2. Desktop GUI      (Tkinter; the processing runs in a background thread)
"""

from __future__ import annotations

import contextlib
import fnmatch
import json
import math
import os
import queue
import re
import subprocess
import sys
import threading
import traceback
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # PDF output uses Agg; the GUI embeds its own Tk canvases
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.ticker import FormatStrFormatter, MultipleLocator  # noqa: E402
from scipy import signal, stats  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

APP_TITLE = "Anisotropy Time-Series Analyzer"
APP_VERSION = "1.0"


# =========================================================================== #
# 1. Processing core
# =========================================================================== #

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

TIME = "Hours"                     # internal name of the time column
TIME_OUT = "Time (h)"              # name of the time column in the Excel output
HOURS_PER_DAY = 24.0
EXCEL_SHEET_NAME_LIMIT = 31
DEFAULT_VALUE_LABEL = "Fluorescence anisotropy (mA)"

REGRESSION_MIN_POINTS = 3           # a line through 2 points has no error estimate
REGRESSION_MIN_PERIOD_HOURS = 12.0  # shorter fitted periods = mis-numbered cycles
ACTOGRAM_MIN_PERIOD_HOURS = 6.0     # range of the adjustable actogram period
ACTOGRAM_MAX_PERIOD_HOURS = 72.0
FIT_MIN_POINTS = 8
FIT_MIN_R2 = 0.1                    # below this a fit is flagged as poor
FIT_DECAY_BOUNDS = (-0.05, 1.0)     # 1/h; negative = growing amplitude


# --------------------------------------------------------------------------- #
# Small utilities
# --------------------------------------------------------------------------- #

def sheet_name(name: str) -> str:
    """Return a sheet name Excel will accept (<= 31 chars, no []:*?/\\)."""
    for bad in "[]:*?/\\":
        name = name.replace(bad, "-")
    return name[:EXCEL_SHEET_NAME_LIMIT]


def safe_filename(text: str) -> str:
    """Make a string usable as part of a file name on Windows/macOS/Linux."""
    return re.sub(r'[\\/:*?"<>|\s]+', "_", str(text).strip()).strip("._")


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def utc_now_string() -> str:
    return utc_now().strftime("%Y-%m-%d %H:%M:%S")


def odd(n: float) -> int:
    n = max(1, int(round(n)))
    return n if n % 2 == 1 else n + 1


def for_excel(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.rename(columns={TIME: TIME_OUT})


def natural_key(text: str):
    """Sort 'A2' before 'A10'."""
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", str(text))]


WELL_RE = re.compile(r"^\s*([A-Za-z]{1,2})\s*0*(\d{1,3})\s*$")


def parse_well(name: str):
    """'D4' -> ('D', 4); None if the name does not look like a plate well."""
    m = WELL_RE.match(str(name))
    return (m.group(1).upper(), int(m.group(2))) if m else None


def match_wells(wells: list[str], pattern: str) -> list[str]:
    """Wells matching a pattern such as 'K*', 'D4, D5 E?' or 'row:K', 'col:5'."""
    tokens = [t for t in re.split(r"[,\s;]+", pattern.strip()) if t]
    chosen = []
    for well in wells:
        parsed = parse_well(well)
        for token in tokens:
            low = token.lower()
            if low.startswith("row:") and parsed:
                if parsed[0] == token[4:].upper():
                    chosen.append(well)
                    break
            elif low.startswith("col:") and parsed:
                if token[4:].isdigit() and parsed[1] == int(token[4:]):
                    chosen.append(well)
                    break
            elif fnmatch.fnmatchcase(well.upper(), token.upper()):
                chosen.append(well)
                break
    return chosen


# --------------------------------------------------------------------------- #
# Reading the data
# --------------------------------------------------------------------------- #

@dataclass
class Recording:
    """One plate-reader time series and the timing facts derived from it."""

    filename: str
    sheet: str
    data: pd.DataFrame            # 'Hours' + one column per well
    temperature: pd.Series        # same index as data (may be all NaN)
    time_header: str = "Time (h)"
    temperature_header: str = "Temperature (\u00b0C)"

    @property
    def wells(self) -> list[str]:
        return [c for c in self.data.columns if c != TIME]

    @property
    def n_points(self) -> int:
        return len(self.data.index)

    @property
    def first_hour(self) -> float:
        return float(self.data[TIME].iloc[0])

    @property
    def last_hour(self) -> float:
        return float(self.data[TIME].iloc[-1])

    @property
    def duration_hours(self) -> float:
        return self.last_hour - self.first_hour

    @property
    def time_interval(self) -> float:
        """Median spacing between samples, in hours (robust to a missed read)."""
        diffs = np.diff(self.data[TIME].to_numpy(dtype=float))
        diffs = diffs[diffs > 0]
        return float(np.median(diffs)) if len(diffs) else 1.0

    def describe(self) -> str:
        return (
            f"{self.n_points} time points, {len(self.wells)} wells, "
            f"{self.first_hour:g}\u2013{self.last_hour:g} h "
            f"({self.duration_hours / HOURS_PER_DAY:.2f} days), "
            f"interval {self.time_interval * 60:.2f} min"
        )

    def subset(self, wells: list[str] | None, start_hour: float, end_hour: float) -> "Recording":
        """Selected wells within [start_hour, end_hour] (end <= start means 'to the end')."""
        wells = [w for w in (wells or self.wells) if w in self.data.columns]
        if not wells:
            raise ValueError("No wells are selected.")
        mask = self.data[TIME] >= start_hour
        if end_hour and end_hour > start_hour:
            mask &= self.data[TIME] <= end_hour
        data = self.data.loc[mask, [TIME] + wells].reset_index(drop=True)
        temperature = self.temperature.loc[mask].reset_index(drop=True)
        return Recording(self.filename, self.sheet, data, temperature,
                         self.time_header, self.temperature_header)


def list_sheets(path) -> list[str]:
    with pd.ExcelFile(path) as book:
        return list(book.sheet_names)


def read_timeseries_excel(path, sheet: str | int | None = None) -> Recording:
    """Read the workbook: col 1 = time (h), col 2 = temperature, col 3.. = wells."""
    frame = pd.read_excel(path, sheet_name=0 if sheet in (None, "") else sheet, header=0)
    frame = frame.dropna(axis=0, how="all")
    if frame.shape[1] < 3:
        raise ValueError(
            "The sheet needs at least 3 columns: time, temperature and one or more wells."
        )
    columns = list(frame.columns)
    hours = pd.to_numeric(frame[columns[0]], errors="coerce")
    keep = hours.notna()
    if keep.sum() < 10:
        raise ValueError(f"Only {int(keep.sum())} numeric time values were found in column 1.")

    data = {TIME: hours[keep].to_numpy(dtype=float)}
    seen: dict[str, int] = {}
    skipped = []
    for column in columns[2:]:
        name = str(column).strip()
        values = pd.to_numeric(frame[column], errors="coerce")[keep].to_numpy(dtype=float)
        if not np.isfinite(values).any():
            skipped.append(name)
            continue
        if name in seen:  # duplicated header -> make it unique
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        data[name] = values
    if len(data) < 2:
        raise ValueError("No numeric well columns were found (column 3 onwards).")

    table = pd.DataFrame(data)
    temperature = pd.Series(
        pd.to_numeric(frame[columns[1]], errors="coerce")[keep].to_numpy(dtype=float),
        name="Temperature",
    )
    order = np.argsort(table[TIME].to_numpy(), kind="stable")  # make sure time increases
    table = table.iloc[order].reset_index(drop=True)
    temperature = temperature.iloc[order].reset_index(drop=True)
    if skipped:
        print(f"Skipped empty/non-numeric columns: {', '.join(skipped)}")

    sheet_label = sheet if isinstance(sheet, str) and sheet else list_sheets(path)[0]
    return Recording(
        filename=str(path),
        sheet=str(sheet_label),
        data=table,
        temperature=temperature,
        time_header=str(columns[0]),
        temperature_header=str(columns[1]),
    )


# --------------------------------------------------------------------------- #
# Pre-processing
# --------------------------------------------------------------------------- #

def remove_outliers(data: pd.DataFrame, wells, window: int, n_mad: float):
    """Hampel-type spike removal: points far from the running median become NaN.

    The scale is the robust SD (1.4826 x MAD) of the residuals from the running
    median over the whole well, so a quiet stretch cannot make every point an
    outlier.  Returns (cleaned data, {well: number removed}).
    """
    cleaned = data.copy()
    removed = {}
    window = odd(window)
    for well in wells:
        series = data[well]
        median = series.rolling(window, center=True, min_periods=1).median()
        residual = series - median
        scale = 1.4826 * np.nanmedian(np.abs(residual - np.nanmedian(residual)))
        if not np.isfinite(scale) or scale <= 0:
            removed[well] = 0
            continue
        bad = residual.abs() > n_mad * scale
        cleaned.loc[bad, well] = np.nan
        removed[well] = int(bad.sum())
    return cleaned, removed


def rolling_mean(data: pd.DataFrame, wells, window: int) -> pd.DataFrame:
    """Centered moving average of the wells; the time column is left alone."""
    smoothed = data[wells].rolling(window=window, center=True, min_periods=1).mean()
    smoothed.insert(0, TIME, data[TIME])
    return smoothed


# --------------------------------------------------------------------------- #
# Detrending
# --------------------------------------------------------------------------- #

DETRENDING_METHODS = ("Sinc Filter", "Moving Average", "Polynomial", "None")


@dataclass
class Trend:
    data: pd.DataFrame
    label: str          # shown in plot legends
    short_label: str    # used in the Excel sheet name
    summary: str        # shown on the 'Note' sheet


def moving_average_trend(data, wells, window_hours: float, interval: float) -> Trend:
    n_points = odd(math.ceil(window_hours / interval))
    print(f"Moving-average trend: {window_hours} h window = {n_points} points")
    trend = rolling_mean(data, wells, n_points)
    trend[wells] = trend[wells].bfill().ffill()
    return Trend(trend, f"{window_hours:g} h moving average ({n_points} points)",
                 f"{n_points}PMA", f"{window_hours:g} h window ({n_points} points)")


def auto_sinc_order(cutoff_hours: float, interval: float, n_points: int) -> int:
    """Taps spanning ~2 x the cutoff period: enough to keep ~24 h rhythms out of a
    48 h trend at any sampling rate, capped by the record length."""
    order = odd(2.0 * cutoff_hours / interval)
    limit = n_points - 1 if (n_points - 1) % 2 == 1 else n_points - 2
    return max(3, min(order, limit))


def sinc_filter_trend(data, wells, cutoff_hours: float, order: int, interval: float) -> Trend:
    n_points = len(data.index)
    nyquist = 0.5 / interval
    norm_cutoff = (1.0 / cutoff_hours) / nyquist
    if not 0 < norm_cutoff < 1:
        raise ValueError(
            f"The sinc cutoff period ({cutoff_hours} h) must be longer than twice the "
            f"sampling interval ({2 * interval:.3f} h)."
        )
    if order <= 0:
        order = auto_sinc_order(cutoff_hours, interval, n_points)
        print(f"Sinc filter order chosen automatically: {order} taps "
              f"({order * interval:.1f} h)")
    order = odd(order)
    limit = n_points - 1 if (n_points - 1) % 2 == 1 else n_points - 2
    if order > limit:
        print(f"Sinc filter order {order} is too long for {n_points} points; using {limit}.")
        order = max(3, limit)
    taps = signal.firwin(order, norm_cutoff, pass_zero="lowpass")
    padlen = min(3 * order, n_points - 1)
    # The ends are extended by point reflection about the mean of the first/last
    # quarter-cutoff of data (not about the single, noisy end point).
    anchor_n = max(3, min(n_points // 4, int(round(cutoff_hours / 4 / interval))))
    print(f"Sinc filter: cutoff {cutoff_hours} h, order {order}, padding {padlen} points")

    trend = data[[TIME] + list(wells)].copy()
    for well in wells:
        if data[well].notna().sum() < 3:
            trend[well] = np.nan
            continue
        y = data[well].interpolate(method="linear", limit_direction="both").to_numpy()
        left = 2 * y[:anchor_n].mean() - y[padlen:0:-1]
        right = 2 * y[-anchor_n:].mean() - y[-2:-padlen - 2:-1]
        extended = np.concatenate([left, y, right])
        try:
            filtered = signal.filtfilt(taps, [1.0], extended, padtype=None)
            trend[well] = filtered[padlen:padlen + n_points]
        except ValueError as error:
            print(f"Warning: could not filter well {well}: {error}")
            trend[well] = np.nan
    return Trend(trend, f"Sinc filter (cutoff {cutoff_hours:g} h, order {order})",
                 f"Sinc C{cutoff_hours:g}h O{order}", f"{cutoff_hours:g} h cutoff, order {order}")


def polynomial_trend(data, wells, degree: int) -> Trend:
    t = data[TIME].to_numpy(dtype=float)
    centre, scale = t.mean(), max(np.ptp(t) / 2, 1e-9)
    x = (t - centre) / scale
    trend = data[[TIME] + list(wells)].copy()
    for well in wells:
        y = data[well].to_numpy(dtype=float)
        ok = np.isfinite(y)
        if ok.sum() <= degree:
            trend[well] = np.nan
            continue
        coefficients = np.polyfit(x[ok], y[ok], degree)
        trend[well] = np.polyval(coefficients, x)
    return Trend(trend, f"Polynomial trend (degree {degree})", f"Poly{degree}",
                 f"Least-squares polynomial of degree {degree}")


def no_trend(data, wells) -> Trend:
    trend = data[[TIME] + list(wells)].copy()
    for well in wells:
        trend[well] = float(np.nanmean(data[well])) if data[well].notna().any() else np.nan
    return Trend(trend, "Mean (no detrending)", "Mean", "None (the mean is subtracted)")


def build_trend(data, wells, method: str, interval: float, s: "Settings") -> Trend:
    if method == "Sinc Filter":
        return sinc_filter_trend(data, wells, s.sinc_cutoff_hours, int(s.sinc_order), interval)
    if method == "Moving Average":
        return moving_average_trend(data, wells, s.ma_window_hours, interval)
    if method == "Polynomial":
        return polynomial_trend(data, wells, int(s.poly_degree))
    if method == "None":
        return no_trend(data, wells)
    raise ValueError(f"Unknown detrending method: {method!r}")


def detrend(data: pd.DataFrame, trend: pd.DataFrame, wells) -> pd.DataFrame:
    detrended = data[wells] - trend[wells]
    detrended.insert(0, TIME, data[TIME])
    return detrended


# --------------------------------------------------------------------------- #
# Peaks and troughs
# --------------------------------------------------------------------------- #

PeakSelection = dict  # well -> (list of peak row labels, list of trough row labels)
RegSelection = dict   # well -> (rows of peaks used, rows of troughs used) in the regression


def find_peaks_and_troughs(series: pd.Series, min_separation_hours: float,
                           interval: float, prominence: float = 0.0):
    valid = series.dropna()
    if len(valid) < 3:
        empty = pd.Index([])
        return empty, empty
    distance = max(1, int(round(min_separation_hours / interval))) if min_separation_hours else 1
    kwargs = {"distance": distance}
    if prominence and prominence > 0:
        kwargs["prominence"] = prominence
    values = valid.to_numpy()
    peaks, _ = signal.find_peaks(values, **kwargs)
    troughs, _ = signal.find_peaks(-values, **kwargs)
    return valid.index[peaks], valid.index[troughs]


def detect_all_peaks(smoothed, wells, s: "Settings", interval: float) -> PeakSelection:
    result = {}
    for well in wells:
        peaks, troughs = find_peaks_and_troughs(
            smoothed[well], s.peak_min_separation_hours, interval, s.peak_min_prominence
        )
        result[well] = ([int(i) for i in peaks], [int(i) for i in troughs])
    return result


def peaks_and_troughs_table(smoothed, wells, selection, auto, reg_selection=None) -> pd.DataFrame:
    rows = []
    for well in wells:
        peaks, troughs = selection[well]
        auto_peaks, auto_troughs = auto[well]
        used = reg_selection.get(well, ([], [])) if reg_selection else (peaks, troughs)
        for kind, indices, automatic, used_rows in (
            ("Peak", peaks, auto_peaks, used[0]),
            ("Trough", troughs, auto_troughs, used[1]),
        ):
            for index in sorted(indices):
                rows.append({
                    "Well": well,
                    "Type": kind,
                    "Time (h)": float(smoothed.loc[index, TIME]),
                    "Smoothed value": float(smoothed.loc[index, well]),
                    "Source": "Auto" if index in automatic else "Manual",
                    "Used in regression": "Yes" if index in set(used_rows) else "No",
                })
    return pd.DataFrame(rows, columns=["Well", "Type", "Time (h)", "Smoothed value",
                                       "Source", "Used in regression"])


def selection_differs(selection: PeakSelection, auto: PeakSelection) -> bool:
    return any(
        sorted(selection[w][0]) != sorted(auto[w][0]) or sorted(selection[w][1]) != sorted(auto[w][1])
        for w in selection
    )


# --------------------------------------------------------------------------- #
# Period / phase regression on the selected peak and trough times
# --------------------------------------------------------------------------- #

REGRESSION_COLUMNS = [
    "Well", "Type", "N Points", "Actogram Period (h)", "Period (h)", "Period SE (h)",
    "R squared", "Residual SD (h)", "Phase (h)", "Phase (degrees)",
    "Fitted Time at Cycle 0 (h)",
]
REGRESSION_POINT_COLUMNS = [
    "Well", "Type", "Cycle", "Time (h)", "Fitted Time (h)", "Residual (h)",
]


def resolve_periods(actogram_period, wells, default: float = HOURS_PER_DAY) -> dict:
    if isinstance(actogram_period, dict):
        return {w: float(actogram_period.get(w) or default) for w in wells}
    value = float(actogram_period) if actogram_period else float(default)
    return {w: value for w in wells}


def default_reg_selection(selection: PeakSelection) -> RegSelection:
    return {w: (list(p), list(t)) for w, (p, t) in selection.items()}


def used_selection(selection_entry, reg_entry) -> tuple[list[int], list[int]]:
    peaks, troughs = selection_entry
    reg_peaks, reg_troughs = reg_entry
    return (sorted(i for i in peaks if i in set(reg_peaks)),
            sorted(i for i in troughs if i in set(reg_troughs)))


def _refine_cycles(gaps: np.ndarray, period: float):
    for _ in range(50):
        steps = np.maximum(1.0, np.rint(gaps / period))
        new_period = float(gaps.sum() / steps.sum())
        if abs(new_period - period) < 1e-9:
            break
        period = new_period
    steps = np.maximum(1.0, np.rint(gaps / period))
    return np.concatenate(([0.0], np.cumsum(steps))), period


def assign_cycle_numbers(times: np.ndarray, period_guess: float) -> np.ndarray:
    """Number the sorted event times 0, 1, 2 ... allowing skipped cycles."""
    gaps = np.diff(times)
    cycles, period = _refine_cycles(gaps, float(period_guess))
    if period < REGRESSION_MIN_PERIOD_HOURS:
        cycles, _ = _refine_cycles(gaps, float(np.median(gaps)))
    return cycles


def regress_period_phase(times, period_guess: float) -> dict | None:
    """Linear regression  time = intercept + period x cycle_number.

    The phase is the fitted event time modulo the period, counted from 0 h.
    """
    t = np.sort(np.asarray(times, dtype=float))
    if len(t) < REGRESSION_MIN_POINTS:
        return None
    cycles = assign_cycle_numbers(t, period_guess)
    if np.ptp(cycles) == 0:
        return None
    fit = stats.linregress(cycles, t)
    period = float(fit.slope)
    if period <= 0:
        return None
    intercept = float(fit.intercept)
    residuals = t - (intercept + period * cycles)
    residual_sd = float(np.sqrt(np.sum(residuals ** 2) / (len(t) - 2))) if len(t) > 2 else np.nan
    phase_hours = intercept % period
    return {
        "n": len(t), "period": period, "period_se": float(fit.stderr),
        "r2": float(fit.rvalue) ** 2, "residual_sd": residual_sd, "intercept": intercept,
        "phase_h": phase_hours, "phase_deg": 360.0 * phase_hours / period,
        "cycles": cycles, "times": t, "fitted": intercept + period * cycles,
    }


def well_regressions(smoothed, peaks, troughs, period_guess: float) -> dict:
    out = {}
    for kind, indices in (("Peak", peaks), ("Trough", troughs)):
        times = smoothed.loc[list(indices), TIME].to_numpy() if len(indices) else []
        out[kind] = regress_period_phase(times, period_guess)
    return out


def compute_all_regressions(smoothed, wells, selection, reg_selection, periods: dict) -> dict:
    result = {}
    for well in wells:
        used = used_selection(selection[well], reg_selection.get(well, selection[well]))
        entry = well_regressions(smoothed, used[0], used[1], periods[well])
        entry["used"] = used
        result[well] = entry
    return result


def regression_summary_text(regression: dict | None, n_peaks: int, n_troughs: int) -> str:
    lines = []
    for kind, name, n in (("Peak", "Peaks", n_peaks), ("Trough", "Troughs", n_troughs)):
        result = regression.get(kind) if regression else None
        if result is None:
            lines.append(f"{name} (n={n}): at least {REGRESSION_MIN_POINTS} selected points "
                         "are needed for the regression")
            continue
        se = result["period_se"]
        se_text = f" \u00b1 {se:.2f} (SE)" if np.isfinite(se) else ""
        lines.append(
            f"{name} (n={n}): period = {result['period']:.2f}{se_text} h,  "
            f"phase = {result['phase_h']:.2f} h ({result['phase_deg']:.1f}\u00b0),  "
            f"R\u00b2 = {result['r2']:.3f}"
        )
    return "\n".join(lines)


def regression_tables(wells, regressions: dict, periods: dict):
    summary_rows, point_rows = [], []
    for well in wells:
        entry = regressions[well]
        for kind, used in (("Peak", entry["used"][0]), ("Trough", entry["used"][1])):
            result = entry[kind]
            row = {"Well": well, "Type": kind, "N Points": len(used),
                   "Actogram Period (h)": periods[well]}
            if result is None:
                row.update({c: np.nan for c in REGRESSION_COLUMNS[4:]})
            else:
                row.update({
                    "Period (h)": result["period"], "Period SE (h)": result["period_se"],
                    "R squared": result["r2"], "Residual SD (h)": result["residual_sd"],
                    "Phase (h)": result["phase_h"], "Phase (degrees)": result["phase_deg"],
                    "Fitted Time at Cycle 0 (h)": result["intercept"],
                })
                for cycle, time, fitted in zip(result["cycles"], result["times"], result["fitted"]):
                    point_rows.append({
                        "Well": well, "Type": kind, "Cycle": int(cycle), "Time (h)": float(time),
                        "Fitted Time (h)": float(fitted), "Residual (h)": float(time - fitted),
                    })
            summary_rows.append(row)
    return (pd.DataFrame(summary_rows, columns=REGRESSION_COLUMNS),
            pd.DataFrame(point_rows, columns=REGRESSION_POINT_COLUMNS))


# --------------------------------------------------------------------------- #
# Saving / loading the hand-edited peaks and troughs
# --------------------------------------------------------------------------- #

PEAKS_JSON_FORMAT = "anisotropy-peaks-v1"


def save_peaks_json(path, selection, hours: pd.Series, source_file: str,
                    reg_selection=None, actogram_period=None) -> None:
    wells = list(selection)
    periods = resolve_periods(actogram_period, wells) if actogram_period is not None else None
    entries = {}
    for well, (peaks, troughs) in selection.items():
        entry = {"peaks": [float(hours.loc[i]) for i in sorted(peaks)],
                 "troughs": [float(hours.loc[i]) for i in sorted(troughs)]}
        if reg_selection is not None and well in reg_selection:
            rp, rt = reg_selection[well]
            entry["regression_peaks"] = [float(hours.loc[i]) for i in sorted(rp)]
            entry["regression_troughs"] = [float(hours.loc[i]) for i in sorted(rt)]
        if periods is not None:
            entry["actogram_period"] = periods[well]
        entries[well] = entry
    payload = {"format": PEAKS_JSON_FORMAT, "source_file": source_file,
               "saved_utc": utc_now_string(), "wells": entries}
    Path(path).write_text(json.dumps(payload, indent=1), encoding="utf-8")


def load_peaks_json(path, hours: pd.Series, wells):
    """Returns (selection, n unmatched times, regression selection, periods)."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("format") != PEAKS_JSON_FORMAT:
        raise ValueError("This is not a peaks/troughs file saved by this application.")
    valid = hours.dropna()
    tolerance = max(1e-3, 0.25 * float(np.median(np.diff(valid.to_numpy())))) if len(valid) > 1 else 1e-3
    skipped = 0

    def to_rows(times) -> list[int]:
        nonlocal skipped
        rows = set()
        for hour in times:
            nearest = (valid - float(hour)).abs().idxmin()
            if abs(float(valid.loc[nearest]) - float(hour)) <= tolerance:
                rows.add(int(nearest))
            else:
                skipped += 1
        return sorted(rows)

    selection, regression, periods = {}, {}, {}
    for well, entry in payload.get("wells", {}).items():
        if well not in wells:
            continue
        selection[well] = (to_rows(entry.get("peaks", [])), to_rows(entry.get("troughs", [])))
        if "regression_peaks" in entry or "regression_troughs" in entry:
            regression[well] = (to_rows(entry.get("regression_peaks", [])),
                                to_rows(entry.get("regression_troughs", [])))
        try:
            value = float(entry.get("actogram_period"))
            if ACTOGRAM_MIN_PERIOD_HOURS <= value <= ACTOGRAM_MAX_PERIOD_HOURS:
                periods[well] = value
        except (TypeError, ValueError):
            pass
    return selection, skipped, (regression or None), (periods or None)


# --------------------------------------------------------------------------- #
# Sine-curve fitting
# --------------------------------------------------------------------------- #

FIT_MODELS = ("Damped sine", "Sine (no damping)")

FIT_COLUMNS = [
    "Well", "Model", "Status", "N Points", "Fit Start (h)", "Fit End (h)",
    "Period (h)", "Period SE (h)",
    "Amplitude", "Amplitude SE", "Relative Amplitude (%)",
    "Peak Phase (h)", "Peak Phase (degrees)", "Phase SE (h)", "Phase (radians)",
    "Decay Rate (1/h)", "Decay Rate SE (1/h)", "Damping Time Constant (h)",
    "Offset", "R squared", "RMSE",
]


def damped_sine(t, amplitude, period, phase, decay, offset):
    """t = time since the start of the fitting window."""
    return amplitude * np.exp(-decay * t) * np.sin(2 * np.pi * t / period + phase) + offset


def plain_sine(t, amplitude, period, phase, offset):
    return amplitude * np.sin(2 * np.pi * t / period + phase) + offset


def _periodogram_guesses(t, y, pmin, pmax, n_guesses: int = 3) -> list[float]:
    """Best periods of a Lomb-Scargle periodogram within [pmin, pmax]."""
    periods = np.linspace(pmin, pmax, 400)
    try:
        power = signal.lombscargle(t, y - y.mean(), 2 * np.pi / periods)
    except Exception:  # noqa: BLE001
        return [float(np.clip(HOURS_PER_DAY, pmin, pmax))]
    peaks, _ = signal.find_peaks(power)
    if len(peaks) == 0:
        peaks = [int(np.argmax(power))]
    best = sorted(peaks, key=lambda i: power[i], reverse=True)[:n_guesses]
    guesses = [float(periods[i]) for i in best]
    mid = float(np.clip(HOURS_PER_DAY, pmin, pmax))
    if all(abs(g - mid) > 1.0 for g in guesses):
        guesses.append(mid)
    return guesses


def _linear_sine(t, y, period):
    """Least-squares amplitude, phase and offset for a fixed period."""
    w = 2 * np.pi / period
    design = np.column_stack([np.sin(w * t), np.cos(w * t), np.ones_like(t)])
    (a, b, c), *_ = np.linalg.lstsq(design, y, rcond=None)
    return float(np.hypot(a, b)), float(np.arctan2(b, a)), float(c)


def fit_sine(times, values, model: str, pmin: float, pmax: float) -> dict:
    """Fit one well.  `times` are absolute hours; returns a dict of results.

    The model uses t' = t - t0 (t0 = first fitted time) so the damping term and
    'Amplitude' refer to the start of the fitting window.  Reported phases are
    converted back to absolute time: 'Peak Phase (h)' is the time of a fitted
    peak modulo the period, counted from 0 h (the same convention as the
    peak-time regression).
    """
    t_abs = np.asarray(times, dtype=float)
    y = np.asarray(values, dtype=float)
    t0 = float(t_abs.min())
    t = t_abs - t0
    damped = model == "Damped sine"
    function = damped_sine if damped else plain_sine
    spread = float(np.nanstd(y)) or 1.0

    best = None
    for period_guess in _periodogram_guesses(t, y, pmin, pmax):
        amplitude, phase, offset = _linear_sine(t, y, period_guess)
        amplitude = max(amplitude, 1e-6 * spread)
        decay_starts = (0.0, 0.03) if damped else (None,)
        for decay in decay_starts:
            if damped:
                p0 = [amplitude * (1 + decay * t.mean()), period_guess, phase, decay, offset]
                lower = [0, pmin, -4 * np.pi, FIT_DECAY_BOUNDS[0], -np.inf]
                upper = [np.inf, pmax, 4 * np.pi, FIT_DECAY_BOUNDS[1], np.inf]
            else:
                p0 = [amplitude, period_guess, phase, offset]
                lower = [0, pmin, -4 * np.pi, -np.inf]
                upper = [np.inf, pmax, 4 * np.pi, np.inf]
            try:
                with np.errstate(over="ignore", invalid="ignore"):
                    params, cov = curve_fit(function, t, y, p0=p0, bounds=(lower, upper),
                                            maxfev=20000)
            except (RuntimeError, ValueError):
                continue
            residual = y - function(t, *params)
            sse = float(np.sum(residual ** 2))
            if np.isfinite(sse) and (best is None or sse < best[0]):
                best = (sse, params, cov)

    if best is None:
        raise RuntimeError("the least-squares fit did not converge")

    sse, params, cov = best
    with np.errstate(invalid="ignore"):
        errors = np.sqrt(np.diag(cov)) if cov is not None else np.full(len(params), np.nan)
    if damped:
        amplitude, period, phase, decay, offset = params
        amplitude_se, period_se, phase_se, decay_se, _ = errors
    else:
        amplitude, period, phase, offset = params
        amplitude_se, period_se, phase_se, _ = errors
        decay, decay_se = 0.0, np.nan

    # absolute phase: sin(2*pi*t/P + phi_abs) with t in hours since 0 h
    phase_abs = (phase - 2 * np.pi * t0 / period) % (2 * np.pi)
    peak_time = (period * (0.25 - phase_abs / (2 * np.pi))) % period
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {
        "params": params, "function": function, "t0": t0,
        "Period (h)": float(period), "Period SE (h)": float(period_se),
        "Amplitude": float(amplitude), "Amplitude SE": float(amplitude_se),
        "Peak Phase (h)": float(peak_time), "Peak Phase (degrees)": float(360 * peak_time / period),
        "Phase SE (h)": float(phase_se * period / (2 * np.pi)) if np.isfinite(phase_se) else np.nan,
        "Phase (radians)": float(phase_abs),
        "Decay Rate (1/h)": float(decay), "Decay Rate SE (1/h)": float(decay_se),
        "Damping Time Constant (h)": float(1.0 / decay) if decay > 1e-9 else np.nan,
        "Offset": float(offset),
        "R squared": 1.0 - sse / ss_tot if ss_tot > 0 else np.nan,
        "RMSE": float(np.sqrt(sse / len(y))),
    }


def fit_window(detrended: pd.DataFrame, start_hour: float, end_hour: float):
    last = float(detrended[TIME].max())
    end = end_hour if end_hour and end_hour > start_hour else last
    return start_hour, end


def fit_one_well(detrended, trend, well, s: "Settings") -> dict:
    """Fit a well's detrended data within the settings' window. Never raises."""
    start, end = fit_window(detrended, s.fit_start_hour, s.fit_end_hour)
    window = detrended[detrended[TIME].between(start, end)]
    values = window[well].dropna()
    times = window.loc[values.index, TIME]
    row = {"Well": well, "Model": s.fit_model, "N Points": len(values),
           "Fit Start (h)": start, "Fit End (h)": end}
    if len(values) < FIT_MIN_POINTS:
        row["Status"] = f"Skipped: only {len(values)} points in the window"
        return {"row": row, "fit": None, "times": times, "values": values}
    try:
        fit = fit_sine(times.to_numpy(), values.to_numpy(), s.fit_model,
                       s.fit_min_period, s.fit_max_period)
    except RuntimeError as error:
        row["Status"] = f"Failed: {error}"
        return {"row": row, "fit": None, "times": times, "values": values}
    level = float(np.nanmean(trend.loc[values.index, well]))
    row.update({k: v for k, v in fit.items() if k in FIT_COLUMNS})
    row["Relative Amplitude (%)"] = 100 * fit["Amplitude"] / abs(level) if level else np.nan
    at_bound = min(abs(fit["Period (h)"] - s.fit_min_period), abs(fit["Period (h)"] - s.fit_max_period))
    if at_bound <= 1e-3:
        row["Status"] = "Period at the search limit"
    elif not fit["R squared"] >= FIT_MIN_R2:
        row["Status"] = f"Poor fit (R\u00b2 < {FIT_MIN_R2:g})"
    else:
        row["Status"] = "OK"
    return {"row": row, "fit": fit, "times": times, "values": values}


def fit_curve(fit: dict, times) -> np.ndarray:
    t = np.asarray(times, dtype=float) - fit["t0"]
    with np.errstate(over="ignore", invalid="ignore"):
        return fit["function"](t, *fit["params"])


def fit_all_wells(detrended, trend, wells, s: "Settings", on_well=None):
    rows, details = [], {}
    for well in wells:
        if on_well is not None:
            on_well(well)
        result = fit_one_well(detrended, trend, well, s)
        if result["fit"] is None:
            print(f"Sine fit, well {well}: {result['row']['Status']}")
        rows.append(result["row"])
        details[well] = result
    table = pd.DataFrame(rows).reindex(columns=FIT_COLUMNS)
    return table, details


# --------------------------------------------------------------------------- #
# Plotting helpers
# --------------------------------------------------------------------------- #

LETTER_PORTRAIT = (8.5, 11.0)   # every PDF page: US Letter, portrait (inches)

OVERVIEW_RC = {"font.size": 5, "axes.titlesize": 6, "axes.labelsize": 5,
               "xtick.labelsize": 4.5, "ytick.labelsize": 4.5, "legend.fontsize": 4}
PAGE_RC = {"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
           "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 6}

OVERVIEW_LAYOUTS = {                 # name: (rows, columns) of graphs per page
    "1 column x 6 rows": (6, 1),
    "2 columns x 6 rows": (6, 2),
    "3 columns x 6 rows": (6, 3),
    "1 column x 8 rows": (8, 1),
    "2 columns x 8 rows": (8, 2),
    "3 columns x 8 rows": (8, 3),
    "4 columns x 8 rows": (8, 4),
}
DEFAULT_OVERVIEW_LAYOUT = "3 columns x 6 rows"
LEGEND_WRAP_CHARS = 20   # maximum characters per line of an outside legend entry


def wrap_label(text: str, width: int = LEGEND_WRAP_CHARS) -> str:
    """Break a legend label into short lines without splitting '48 h' or 'T = 24 h'."""
    import textwrap

    nbsp = "\u00a0"
    text = re.sub(r"(\d) (h|min|points?)\b", rf"\1{nbsp}\2", str(text))
    text = re.sub(r"\b(order|degree|cutoff) (\d)", rf"\1{nbsp}\2", text)
    text = text.replace(" = ", f"{nbsp}={nbsp}")
    lines = textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False)
    return "\n".join(lines).replace(nbsp, " ")


def outside_legend(ax, fontsize: float = 6.5, width: int = LEGEND_WRAP_CHARS, **kwargs):
    """Legend to the right of the axes (top-aligned), with wrapped labels."""
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    legend = ax.legend(handles, [wrap_label(label, width) for label in labels],
                       loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0,
                       fontsize=fontsize, handlelength=1.6, labelspacing=0.6, **kwargs)
    for handle in getattr(legend, "legend_handles", getattr(legend, "legendHandles", [])):
        # data-point symbols are drawn tiny in the plot; make them visible in the legend
        if hasattr(handle, "get_sizes") and len(handle.get_sizes()) and max(handle.get_sizes()) < 12:
            handle.set_sizes([12])
        elif hasattr(handle, "get_markersize") and handle.get_markersize() < 3.5:
            handle.set_markersize(3.5)
    return legend


MAJOR_TICK_HOURS = {"Every 12 hours": 12, "Every 24 hours": 24, "Every 48 hours": 48}
MINOR_TICK_HOURS = {"No minor ticks": 0, "Every 2 hours": 2, "Every 4 hours": 4,
                    "Every 6 hours": 6, "Every 12 hours": 12}


@dataclass
class AxisStyle:
    x_min: float
    x_max: float
    major: int
    minor: float

    def apply(self, ax) -> None:
        ax.set_xlim(self.x_min, self.x_max)
        ax.xaxis.set_major_locator(MultipleLocator(self.major))
        if self.minor and self.minor < self.major:
            ax.xaxis.set_minor_locator(MultipleLocator(self.minor))
        ax.grid(True, linewidth=0.5, color="lightgray", linestyle="--")


def make_axis_style(hours: pd.Series, major: int, minor: float) -> AxisStyle:
    lo, hi = float(hours.min()), float(hours.max())
    x_min = math.floor(lo / major) * major if lo >= 0 else lo
    x_max = math.ceil(hi / major) * major
    if x_max - hi < 0.02 * (hi - lo):
        x_max += major / 2
    return AxisStyle(x_min, x_max, major, minor)


def overview_pages(wells: list[str], layout: str):
    """[(n_rows, n_cols, [(row, col, well), ...]), ...] - one entry per PDF page.

    The wells are placed in file order, left to right and then top to bottom.
    """
    n_rows, n_cols = OVERVIEW_LAYOUTS.get(layout, OVERVIEW_LAYOUTS[DEFAULT_OVERVIEW_LAYOUT])
    per_page = n_rows * n_cols
    pages = []
    for first in range(0, len(wells), per_page):
        chunk = wells[first:first + per_page]
        pages.append((n_rows, n_cols, [(i // n_cols, i % n_cols, w) for i, w in enumerate(chunk)]))
    return pages


def plot_overview(pdf, *, title, wells, layout, scatter, line, line_label, y_label, axis) -> None:
    pages = overview_pages(wells, layout)
    for page_number, (n_rows, n_cols, cells) in enumerate(pages, start=1):
        # wider graphs (1-2 columns) get larger text and points
        scale = {1: 1.5, 2: 1.25}.get(n_cols, 1.0)
        rc = {key: value * scale for key, value in OVERVIEW_RC.items()}
        with plt.rc_context(rc):
            figure = plt.figure(figsize=LETTER_PORTRAIT)
            suffix = f"  (page {page_number}/{len(pages)})" if len(pages) > 1 else ""
            figure.suptitle(title + suffix, fontsize=11)
            grid = figure.add_gridspec(n_rows, n_cols, hspace=0.45, wspace=0.28,
                                       left=0.07, right=0.98, top=0.93, bottom=0.07)
            for r, c, well in cells:
                ax = figure.add_subplot(grid[r, c])
                ax.scatter(scatter[TIME], scatter[well], s=0.3 * scale ** 2, c="tab:blue",
                           linewidths=0)
                ax.plot(line[TIME], line[well], "-r", linewidth=0.6 * scale, label=line_label)
                axis.apply(ax)
                ax.set_title(well, pad=2)
            figure.text(0.52, 0.015, "Time (h)", ha="center", fontsize=9)
            figure.text(0.012, 0.5, y_label, va="center", rotation="vertical", fontsize=9)
            figure.text(0.98, 0.015, f"red line: {line_label}", ha="right", fontsize=6,
                        color="red")
            pdf.savefig(figure)
            plt.close(figure)


def plot_temperature_page(pdf, recording: Recording, title: str, axis: AxisStyle) -> None:
    """Temperature time course and its distribution (US Letter, portrait)."""
    temp = recording.temperature
    if not temp.notna().any():
        return
    values = temp.dropna().to_numpy(dtype=float)
    with plt.rc_context(PAGE_RC):
        figure = plt.figure(figsize=LETTER_PORTRAIT)
        figure.suptitle(f"{title} \u2013 temperature", fontsize=12)
        grid = figure.add_gridspec(2, 1, height_ratios=[1.25, 1.0], hspace=0.35,
                                   left=0.12, right=0.95, top=0.92, bottom=0.08)
        ax = figure.add_subplot(grid[0])
        ax.plot(recording.data[TIME], temp, "-", color="tab:orange", linewidth=0.8)
        axis.apply(ax)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Temperature (\u00b0C)")
        ax.set_title(
            f"Mean {values.mean():.2f} \u00b0C,  SD {values.std():.2f} \u00b0C,  "
            f"range {values.min():.1f}\u2013{values.max():.1f} \u00b0C  (n = {len(values)})",
            fontsize=10)

        ax_hist = figure.add_subplot(grid[1])
        levels = np.unique(values)
        if 1 < len(levels) <= 60:
            # the reader reports temperature in fixed steps (e.g. 0.1 degC): one bar per step
            step = float(np.diff(levels).min())
            edges = np.arange(levels[0] - step / 2, levels[-1] + step, step)
        else:
            edges = 40 if len(levels) > 1 else 1
        ax_hist.hist(values, bins=edges, color="tab:orange", alpha=0.8, edgecolor="white")
        ax_hist.yaxis.get_major_formatter().set_useOffset(False)
        ax_hist.axvline(values.mean(), color="black", linestyle="--", linewidth=1,
                        label=f"mean {values.mean():.2f} \u00b0C")
        ax_hist.set_xlabel("Temperature (\u00b0C)")
        ax_hist.set_ylabel("Number of time points")
        ax_hist.set_title("Distribution of the recorded temperatures", fontsize=10)
        ax_hist.legend(fontsize=8)
        ax_hist.grid(True, linestyle="--", alpha=0.5)
        pdf.savefig(figure)
        plt.close(figure)


def actogram_points(smoothed: pd.DataFrame, indices, day_length: float) -> list:
    """(x, row, hour, row label) of every drawn copy of each event (double plot)."""
    points = []
    if len(indices) == 0:
        return points
    hours = smoothed.loc[list(indices), TIME]
    for index, hour in zip(hours.index, hours.to_numpy()):
        day = hour // day_length
        time_in_day = hour - day * day_length
        points.append((time_in_day, day, hour, int(index)))
        if day > 0:
            points.append((time_in_day + day_length, day - 1, hour, int(index)))
    return points


def _regression_line(result: dict, day_length: float):
    cycles = np.arange(0, int(result["cycles"].max()) + 1)
    times = result["intercept"] + result["period"] * cycles
    day = np.floor(times / day_length)
    x = times - day * day_length
    row = day.copy()
    for i in range(1, len(x)):
        while x[i] - x[i - 1] > day_length / 2:
            x[i] -= day_length
            row[i] += 1
        while x[i] - x[i - 1] < -day_length / 2:
            x[i] += day_length
            row[i] -= 1
    return x, row


def draw_actogram(ax, smoothed, peaks, troughs, day_length: float, first_hour: float,
                  last_hour: float, show_labels: bool, selected=None, regression=None,
                  summary: str = "", legend_fontsize: float = 6.5) -> None:
    """Double-plotted actogram of peak/trough times (ringed = used in regression)."""
    selected_sets = (set(selected[0]), set(selected[1])) if selected is not None else None
    ring_labelled = False
    for kind_number, (indices, color, dark, label) in enumerate(
        ((peaks, "red", "darkred", "Peaks"), (troughs, "blue", "navy", "Troughs"))
    ):
        points = actogram_points(smoothed, indices, day_length)
        if points:
            xs, ys, original_hours, rows = zip(*points)
            ax.scatter(xs, ys, marker="o", s=20, color=color, label=label, zorder=3)
            if show_labels:
                for x, y, hour in zip(xs, ys, original_hours):
                    ax.text(x + 0.02 * day_length, y, f"{hour:.1f}", fontsize=5, color=color)
            if selected_sets is not None:
                ringed = [(x, y) for x, y, _h, r in points if r in selected_sets[kind_number]]
                if ringed:
                    rx, ry = zip(*ringed)
                    ax.scatter(rx, ry, marker="o", s=85, facecolors="none", edgecolors="black",
                               linewidths=1.1, zorder=4,
                               label=None if ring_labelled else "Used in regression")
                    ring_labelled = True
        result = regression.get("Peak" if kind_number == 0 else "Trough") if regression else None
        if result is not None:
            line_x, line_row = _regression_line(result, day_length)
            for k, shift in enumerate((0, 1, -1)):
                ax.plot(line_x + shift * day_length, line_row - shift, linestyle="--",
                        linewidth=1.0, color=dark, alpha=0.9, zorder=2,
                        label=(f"{label[:-1]} regression (T = {result['period']:.2f} h)"
                               if k == 0 else None))

    pad = 0.085 if show_labels else 0.0
    ax.set_xlim(-pad * day_length, (2 + pad) * day_length)
    ax.set_xticks([f * day_length for f in np.arange(0, 2.25, 0.25)])
    decimals = 0 if (day_length / 4).is_integer() else (1 if (day_length * 10 / 4).is_integer() else 2)
    ax.xaxis.set_major_formatter(FormatStrFormatter(f"%.{decimals}f"))
    ax.set_title(f"Actogram of peaks and troughs (double-plotted, T = {day_length:g} h)",
                 fontsize=10)
    ax.set_ylabel(f"Cycle (T = {day_length:g} h)")
    ax.set_xlabel("Time within cycle (h)")
    first_row = max(first_hour, 0) // day_length
    n_rows = max(last_hour // day_length, first_row + 1)
    ax.set_ylim(first_row - 0.6, n_rows + 0.6)
    ax.yaxis.set_major_locator(MultipleLocator(base=1))
    ax.grid(True, linestyle="--", alpha=0.7)
    outside_legend(ax, fontsize=legend_fontsize)
    ax.invert_yaxis()
    if summary:
        ax.text(0.5, -0.24, summary, transform=ax.transAxes, ha="center", va="top",
                fontsize=8.5, linespacing=1.5,
                bbox={"boxstyle": "round,pad=0.4", "facecolor": "#f4f4f4", "edgecolor": "gray"})


def fit_summary_text(row: dict) -> str:
    if row.get("Status", "").startswith(("Skipped", "Failed")):
        return f"Sine fit: {row['Status']}"
    text = (f"Sine fit ({row['Model']}): period = {row['Period (h)']:.2f} h, "
            f"amplitude = {row['Amplitude']:.3g}, peak phase = {row['Peak Phase (h)']:.2f} h "
            f"({row['Peak Phase (degrees)']:.1f}\u00b0), R\u00b2 = {row['R squared']:.3f}")
    if row["Model"] == "Damped sine":
        text += f", decay = {row['Decay Rate (1/h)']:.4f} /h"
    if row.get("Status") != "OK":
        text += f"  [{row['Status']}]"
    return text


def plot_well_page(pdf, *, well, title, value_label, raw, trend, detrended, smoothed, axis,
                   smooth_label, peaks_troughs, reg_selection, regression, day_length,
                   first_hour, last_hour, label_peaks, label_actogram) -> None:
    """Raw + trend, detrended + peaks/troughs, actogram with regression."""
    with plt.rc_context(PAGE_RC):
        figure = plt.figure(figsize=LETTER_PORTRAIT)
        figure.suptitle(f"{title}   Well {well}", fontsize=12)
        grid = figure.add_gridspec(3, 1, height_ratios=[1, 1, 1.25], hspace=0.38,
                                   left=0.10, right=0.80, top=0.94, bottom=0.12)

        ax_raw = figure.add_subplot(grid[0])
        ax_raw.scatter(raw[TIME], raw[well], s=3.0, c="tab:purple", linewidths=0, label="Raw data")
        ax_raw.plot(trend.data[TIME], trend.data[well], "-r", linewidth=1.0,
                    label=f"Trend: {trend.label}")
        axis.apply(ax_raw)
        ax_raw.set_ylabel(value_label)
        ax_raw.set_xlabel("Time (h)")
        outside_legend(ax_raw)

        ax_det = figure.add_subplot(grid[1])
        ax_det.scatter(detrended[TIME], detrended[well], s=3.0, c="tab:purple", linewidths=0,
                       label="Detrended")
        ax_det.plot(smoothed[TIME], smoothed[well], "-b", linewidth=1.0, label=smooth_label)
        peaks, troughs = peaks_troughs
        low, high = ax_det.get_ylim()
        span = high - low
        for indices, color, label, sign, va in ((peaks, "red", "Peaks", 1, "bottom"),
                                                (troughs, "blue", "Troughs", -1, "top")):
            if not indices:
                continue
            hours = smoothed.loc[indices, TIME]
            values = smoothed.loc[indices, well]
            ax_det.scatter(hours, values, marker="o", s=30, color=color, label=label, zorder=4)
            if label_peaks:
                for hour, value in zip(hours, values):
                    ax_det.text(hour, value + sign * span * 0.05, f"{hour:.1f}", fontsize=5,
                                color=color, ha="center", va=va)
        axis.apply(ax_det)
        ax_det.set_ylabel(f"Detrended {value_label.split('(')[0].strip().lower()}")
        ax_det.set_xlabel("Time (h)")
        outside_legend(ax_det)

        ax_act = figure.add_subplot(grid[2])
        summary = regression_summary_text(regression, len(reg_selection[0]), len(reg_selection[1]))
        draw_actogram(ax_act, smoothed, peaks, troughs, day_length, first_hour, last_hour,
                      label_actogram, selected=reg_selection, regression=regression,
                      summary=summary)
        pdf.savefig(figure)
        plt.close(figure)


def plot_fit_page(pdf, *, well, title, value_label, raw, trend, detrended, detail, axis) -> None:
    fit, row = detail["fit"], detail["row"]
    times, values = detail["times"], detail["values"]
    with plt.rc_context(PAGE_RC):
        figure, (ax_raw, ax_det) = plt.subplots(2, 1, figsize=LETTER_PORTRAIT, sharex=True)
        figure.suptitle(f"{title} \u2013 {row['Model']} fit, well {well}", fontsize=12)
        curve = fit_curve(fit, times)
        ax_raw.plot(raw[TIME], raw[well], "o", markersize=1.8, color="gray", label="Raw data")
        ax_raw.plot(trend.data[TIME], trend.data[well], "k--", linewidth=1.0, label="Trend")
        ax_raw.plot(times, curve + trend.data.loc[times.index, well], "r-", linewidth=1.4,
                    label="Fit + trend")
        ax_raw.set_ylabel(value_label)
        ax_raw.set_title(f"P = {row['Period (h)']:.2f} h,  A = {row['Amplitude']:.3g},  "
                         f"peak phase = {row['Peak Phase (h)']:.2f} h", fontsize=10)
        outside_legend(ax_raw)
        axis.apply(ax_raw)

        ax_det.plot(detrended[TIME], detrended[well], "o", markersize=1.8, color="lightgray",
                    label="Detrended (all)")
        ax_det.plot(times, values, "o", markersize=1.8, color="tab:blue",
                    label=f"Fitted range {row['Fit Start (h)']:.1f}\u2013{row['Fit End (h)']:.1f} h")
        ax_det.plot(times, curve, "r-", linewidth=1.4, label="Fit")
        ax_det.set_xlabel("Time (h)")
        ax_det.set_ylabel("Detrended value")
        extra = (f",  decay = {row['Decay Rate (1/h)']:.4f} /h"
                 if row["Model"] == "Damped sine" else "")
        ax_det.set_title(f"R\u00b2 = {row['R squared']:.3f},  RMSE = {row['RMSE']:.3g}{extra}",
                         fontsize=10)
        outside_legend(ax_det)
        axis.apply(ax_det)
        figure.tight_layout(rect=(0, 0.02, 1, 0.96))
        pdf.savefig(figure)
        plt.close(figure)


PLATE_ROWS = "ABCDEFGHIJKLMNOP"      # 384-well plate: 16 rows x 24 columns
PLATE_COLUMNS = 24
PLATE_DATA_COLOR = "#b9ecb9"        # light green: wells with data
PLATE_EMPTY_COLOR = "white"


def plate_position(well: str):
    """(row index 0-15, column index 0-23) of a 384-well plate, or None."""
    parsed = parse_well(well)
    if parsed is None:
        return None
    letter, number = parsed
    if len(letter) != 1 or letter not in PLATE_ROWS or not 1 <= number <= PLATE_COLUMNS:
        return None
    return PLATE_ROWS.index(letter), number - 1


def plot_plate_map(pdf, data_wells, analysed_wells, title: str) -> None:
    """384-well plate map (US Letter, portrait); wells with data are light green.

    Dimensions follow the ANSI/SLAS microplate standard (127.76 x 85.48 mm,
    A1 centre 12.13 mm from the left and 8.99 mm from the top, 4.5 mm pitch).
    If only some of the wells were analysed, those get a dark green outline.
    """
    from matplotlib.patches import FancyBboxPatch

    width, height = 127.76, 85.48
    x_a1, y_a1, pitch, well_size = 12.13, 8.99, 4.5, 3.7
    positions = {w: plate_position(w) for w in data_wells}
    on_plate = {w: p for w, p in positions.items() if p is not None}
    off_plate = [w for w, p in positions.items() if p is None]
    analysed = set(analysed_wells)
    subset = analysed != set(data_wells)

    with plt.rc_context(PAGE_RC):
        figure = plt.figure(figsize=LETTER_PORTRAIT)
        figure.suptitle(f"{title} \u2013 plate map (384-well plate)", fontsize=12)
        ax = figure.add_axes([0.05, 0.47, 0.90, 0.42])
        ax.set_xlim(-4, width + 1)
        ax.set_ylim(height + 1, -6)          # row A at the top
        ax.set_aspect("equal")
        ax.axis("off")
        ax.add_patch(FancyBboxPatch((0, 0), width, height, boxstyle="round,pad=0,rounding_size=3",
                                    facecolor="#f4f4f4", edgecolor="#555555", linewidth=1.2))
        for c in range(PLATE_COLUMNS):
            ax.text(x_a1 + c * pitch, y_a1 - 4.2, str(c + 1), ha="center", va="center", fontsize=6)
        for r, letter in enumerate(PLATE_ROWS):
            ax.text(x_a1 - 4.4, y_a1 + r * pitch, letter, ha="center", va="center", fontsize=6)

        filled = {p: w for w, p in on_plate.items()}
        for r in range(len(PLATE_ROWS)):
            for c in range(PLATE_COLUMNS):
                x = x_a1 + c * pitch - well_size / 2
                y = y_a1 + r * pitch - well_size / 2
                well = filled.get((r, c))
                has_data = well is not None
                edge, lw = "#9a9a9a", 0.5
                if has_data and subset and well in analysed:
                    edge, lw = "darkgreen", 1.3
                ax.add_patch(FancyBboxPatch(
                    (x, y), well_size, well_size, boxstyle="round,pad=0,rounding_size=0.5",
                    facecolor=PLATE_DATA_COLOR if has_data else PLATE_EMPTY_COLOR,
                    edgecolor=edge, linewidth=lw))
                if has_data:
                    ax.text(x + well_size / 2, y + well_size / 2, well, ha="center", va="center",
                            fontsize=4.2, color="#1f5f1f")

        # legend and a short description below the plate
        ax_key = figure.add_axes([0.08, 0.05, 0.84, 0.38])
        ax_key.axis("off")
        ax_key.set_xlim(0, 1)
        ax_key.set_ylim(0, 1)

        def key(y, face, edge, lw, text):
            ax_key.add_patch(FancyBboxPatch((0.0, y - 0.022), 0.03, 0.044,
                                            boxstyle="round,pad=0,rounding_size=0.004",
                                            facecolor=face, edgecolor=edge, linewidth=lw,
                                            transform=ax_key.transAxes))
            ax_key.text(0.045, y, text, va="center", fontsize=9)

        key(0.95, PLATE_DATA_COLOR, "#9a9a9a", 0.5, f"Well with data ({len(on_plate)})")
        y = 0.87
        if subset:
            key(y, PLATE_DATA_COLOR, "darkgreen", 1.3,
                f"Analysed in this run ({len(analysed & set(on_plate))})")
            y -= 0.08
        key(y, PLATE_EMPTY_COLOR, "#9a9a9a", 0.5,
            f"No data ({len(PLATE_ROWS) * PLATE_COLUMNS - len(on_plate)})")
        y -= 0.10
        if on_plate:
            rows = sorted({p[0] for p in on_plate.values()})
            cols = sorted({p[1] for p in on_plate.values()})
            ax_key.text(0.0, y, f"Rows {PLATE_ROWS[rows[0]]}\u2013{PLATE_ROWS[rows[-1]]}, "
                                f"columns {cols[0] + 1}\u2013{cols[-1] + 1}", fontsize=9, va="top")
            y -= 0.07
        if off_plate:
            names = ", ".join(off_plate[:40]) + (" ..." if len(off_plate) > 40 else "")
            ax_key.text(0.0, y, f"Not shown (not a 384-well position): {names}",
                        fontsize=8, va="top", wrap=True)
        pdf.savefig(figure)
        plt.close(figure)


def plot_summary_pages(pdf, summary: pd.DataFrame, title: str, s: "Settings",
                       per_page: int = 48) -> None:
    """Period, phase and amplitude of every well (regression and sine fit).

    Filled black circles: sine fits with status OK; open grey circles: other fits.
    """
    for first in range(0, len(summary), per_page):
        part = summary.iloc[first:first + per_page].reset_index(drop=True)
        x = np.arange(len(part))
        ok = (part["Fit Status"] == "OK").to_numpy()
        with plt.rc_context(PAGE_RC):
            figure, axes = plt.subplots(3, 1, figsize=LETTER_PORTRAIT, sharex=True)
            suffix = f" (wells {first + 1}\u2013{first + len(part)})" if len(summary) > per_page else ""
            figure.suptitle(f"{title} \u2013 summary{suffix}", fontsize=12)

            ax = axes[0]
            period = part["Fit Period (h)"].to_numpy(dtype=float)
            err = part["Fit Period SE (h)"].to_numpy(dtype=float)
            ax.errorbar(x[ok], period[ok], yerr=err[ok], fmt="o", color="black", ms=4,
                        capsize=2, label="Sine fit")
            ax.plot(x[~ok], period[~ok], "o", mfc="none", mec="gray", ms=4,
                    label="Sine fit (not OK)")
            ax.plot(x - 0.18, part["Peak Regression Period (h)"], "^", color="red", ms=4,
                    label="Peak regression")
            ax.plot(x + 0.18, part["Trough Regression Period (h)"], "v", color="blue", ms=4,
                    label="Trough regression")
            ax.set_ylim(s.fit_min_period * 0.9, s.fit_max_period * 1.05)
            ax.set_ylabel("Period (h)")
            outside_legend(ax)

            ax = axes[1]
            phase = part["Fit Peak Phase (h)"].to_numpy(dtype=float)
            ax.plot(x[ok], phase[ok], "o", color="black", ms=4, label="Sine fit")
            ax.plot(x[~ok], phase[~ok], "o", mfc="none", mec="gray", ms=4,
                    label="Sine fit (not OK)")
            ax.plot(x - 0.18, part["Peak Regression Phase (h)"], "^", color="red", ms=4,
                    label="Peak regression")
            ax.set_ylim(bottom=0)
            ax.set_ylabel("Peak phase (h, from 0 h)")
            outside_legend(ax)

            ax = axes[2]
            amp = part["Fit Amplitude"].to_numpy(dtype=float)
            ax.bar(x, np.where(ok, amp, np.nan), color="tab:purple", alpha=0.75, label="OK")
            ax.bar(x, np.where(ok, np.nan, amp), color="lightgray", label="not OK")
            ax.set_ylabel("Fit amplitude")
            outside_legend(ax)
            ax.set_xticks(x)
            ax.set_xticklabels(part["Well"], rotation=90, fontsize=6.5)
            ax.set_xlim(-0.8, max(len(part), 12) - 0.2)
            for a in axes:
                a.grid(True, linestyle="--", alpha=0.5)
            figure.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(figure)
            plt.close(figure)


# --------------------------------------------------------------------------- #
# Settings, analysis and output
# --------------------------------------------------------------------------- #

@dataclass
class Settings:
    experiment_name: str = "experiment"
    experiment_title: str = ""
    value_label: str = DEFAULT_VALUE_LABEL
    sheet: str = ""
    wells: list = field(default_factory=list)       # empty = all wells
    start_hour: float = 0.0
    end_hour: float = 0.0                            # 0 = to the end
    remove_outliers: bool = False
    outlier_window: int = 7
    outlier_threshold: float = 5.0
    detrending_method: str = "Sinc Filter"
    sinc_cutoff_hours: float = 48.0
    sinc_order: int = 0                              # 0 = automatic
    ma_window_hours: float = 24.0
    poly_degree: int = 3
    smoothing_points: int = 9
    peak_min_separation_hours: float = 12.0
    peak_min_prominence: float = 0.0
    overview_layout: str = DEFAULT_OVERVIEW_LAYOUT
    major_ticks: str = "Every 24 hours"
    minor_ticks: str = "Every 12 hours"
    label_detrended_plot: bool = True
    label_actogram: bool = False
    actogram_period: float = 24.0
    well_pages: bool = True
    fit_pages: bool = True
    well_sheets: bool = True
    fit_model: str = "Damped sine"
    fit_start_hour: float = 12.0
    fit_end_hour: float = 0.0                        # 0 = to the end
    fit_min_period: float = 12.0
    fit_max_period: float = 48.0


def expected_outputs(output_dir, settings: Settings) -> list[Path]:
    name = safe_filename(settings.experiment_name) or "experiment"
    d = Path(output_dir)
    return [d / f"{name}_results.xlsx", d / f"{name}_data.xlsx",
            d / f"{name}_plots.pdf", d / f"{name}_peaks.json"]


@dataclass
class Analysis:
    input_path: Path
    settings: Settings
    source: Recording        # whole file
    recording: Recording     # selected wells / time range (after outlier removal)
    wells: list
    time_interval: float
    outliers: dict
    trend: Trend
    detrended: pd.DataFrame
    smoothed: pd.DataFrame           # smoothed detrended data (peak detection)
    smoothed_raw: pd.DataFrame
    auto_peaks: PeakSelection

    @property
    def raw(self) -> pd.DataFrame:
        return self.recording.data

    @property
    def smooth_label(self) -> str:
        n = int(self.settings.smoothing_points)
        return f"{n}-point moving average ({(n - 1) * self.time_interval:.2g} h)"


def _reporter(progress):
    def report(fraction: float, text: str) -> None:
        if progress is not None:
            progress(min(max(fraction, 0.0), 1.0), text)
    return report


def analyze(input_path, settings: Settings, progress=None, source: Recording | None = None) -> Analysis:
    """Read, clean, detrend, smooth and detect peaks/troughs (no files written)."""
    report = _reporter(progress)
    s = settings
    input_path = Path(input_path)
    print(f"Started at {utc_now_string()} (UTC)")
    report(0.0, "Reading the data...")
    if source is None:
        print(f"Reading {input_path.name} ...")
        source = read_timeseries_excel(input_path, s.sheet or None)
    print(f"File: {source.describe()}")
    recording = source.subset(s.wells, s.start_hour, s.end_hour)
    wells = recording.wells
    if recording.n_points < 10:
        raise ValueError(f"Only {recording.n_points} time points are inside the analysis range.")
    interval = recording.time_interval
    print(f"Analysis: {recording.describe()}")

    outliers = {}
    if s.remove_outliers:
        report(0.15, "Removing outliers...")
        cleaned, outliers = remove_outliers(recording.data, wells, s.outlier_window,
                                            s.outlier_threshold)
        recording = replace(recording, data=cleaned)
        print(f"Outlier removal: {sum(outliers.values())} points set to blank "
              f"({', '.join(f'{w}:{n}' for w, n in outliers.items() if n) or 'none'})")

    report(0.3, "Detrending...")
    print(f"Detrending with: {s.detrending_method}")
    trend = build_trend(recording.data, wells, s.detrending_method, interval, s)
    detrended = detrend(recording.data, trend.data, wells)
    window = odd(s.smoothing_points)
    smoothed = rolling_mean(detrended, wells, window)
    smoothed_raw = rolling_mean(recording.data, wells, window)

    report(0.8, "Detecting peaks and troughs...")
    auto_peaks = detect_all_peaks(smoothed, wells, s, interval)
    n_p = sum(len(p) for p, _ in auto_peaks.values())
    n_t = sum(len(t) for _, t in auto_peaks.values())
    print(f"Detected {n_p} peaks and {n_t} troughs in {len(wells)} wells.")
    report(1.0, "Analysis finished")
    return Analysis(input_path, s, source, recording, wells, interval, outliers, trend,
                    detrended, smoothed, smoothed_raw, auto_peaks)


def build_summary(wells, regressions, fit_table: pd.DataFrame, recording: Recording) -> pd.DataFrame:
    fits = fit_table.set_index("Well")
    rows = []
    for well in wells:
        entry = regressions[well]
        row = {"Well": well, "Mean Raw Value": float(np.nanmean(recording.data[well]))}
        for kind in ("Peak", "Trough"):
            result = entry[kind]
            n_used = len(entry["used"][0 if kind == "Peak" else 1])
            row[f"{kind} Regression N"] = n_used
            row[f"{kind} Regression Period (h)"] = result["period"] if result else np.nan
            row[f"{kind} Regression Period SE (h)"] = result["period_se"] if result else np.nan
            row[f"{kind} Regression Phase (h)"] = result["phase_h"] if result else np.nan
            row[f"{kind} Regression R squared"] = result["r2"] if result else np.nan
        fit = fits.loc[well]
        for key in ("Period (h)", "Period SE (h)", "Amplitude", "Relative Amplitude (%)",
                    "Peak Phase (h)", "Peak Phase (degrees)", "Decay Rate (1/h)", "R squared"):
            row[f"Fit {key}"] = fit.get(key, np.nan)
        row["Fit Status"] = fit.get("Status", "")
        rows.append(row)
    return pd.DataFrame(rows)


def build_note_sheet(rows) -> pd.DataFrame:
    return pd.DataFrame({"Item": [r[0] for r in rows], "Value": [r[1] for r in rows]})


def autosize(writer, sheet: str, frame: pd.DataFrame) -> None:
    """Reasonable column widths (openpyxl engine)."""
    try:
        ws = writer.sheets[sheet]
    except KeyError:
        return
    from openpyxl.utils import get_column_letter
    for i, column in enumerate(frame.columns, start=1):
        width = max(len(str(column)), 8)
        ws.column_dimensions[get_column_letter(i)].width = min(width + 2, 45)
    ws.freeze_panes = "B2"


def write_outputs(analysis: Analysis, output_dir, peaks=None, progress=None,
                  reg_selection=None, actogram_period=None) -> list[Path]:
    """Excel workbooks, PDF and the peaks JSON.  Returns the files written."""
    report = _reporter(progress)
    s = analysis.settings
    wells = analysis.wells
    raw, trend = analysis.raw, analysis.trend
    detrended, smoothed = analysis.detrended, analysis.smoothed
    peaks = peaks or analysis.auto_peaks
    reg_selection = reg_selection or default_reg_selection(peaks)
    edited = selection_differs(peaks, analysis.auto_peaks)
    periods = resolve_periods(actogram_period if actogram_period is not None else s.actogram_period,
                              wells, s.actogram_period)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_xlsx, data_xlsx, plots_pdf, peaks_json = expected_outputs(output_dir, s)
    start_time = utc_now()
    title = s.experiment_name + (f" \u2013 {s.experiment_title}" if s.experiment_title else "")

    report(0.02, "Period/phase regression...")
    regressions = compute_all_regressions(smoothed, wells, peaks, reg_selection, periods)
    reg_table, reg_points = regression_tables(wells, regressions, periods)

    report(0.05, "Sine fitting...")
    print(f"\nSine-curve fitting ({s.fit_model})...")
    counter = {"n": 0}

    def on_well(well):
        counter["n"] += 1
        report(0.05 + 0.20 * counter["n"] / len(wells), f"Fitting well {well}...")

    fit_table, fit_details = fit_all_wells(detrended, trend.data, wells, s, on_well)
    summary = build_summary(wells, regressions, fit_table, analysis.recording)
    table = peaks_and_troughs_table(smoothed, wells, peaks, analysis.auto_peaks,
                                    {w: regressions[w]["used"] for w in wells})
    fit_start, fit_end = fit_window(detrended, s.fit_start_hour, s.fit_end_hour)

    try:
        n_sheets = len(list_sheets(analysis.input_path))
    except Exception:  # noqa: BLE001 - only decides whether the sheet name is noted
        n_sheets = 1
    note = build_note_sheet([
        ("Application", f"{APP_TITLE} {APP_VERSION}"),
        ("Experiment", s.experiment_name),
        ("Title", s.experiment_title),
        ("Input file", analysis.input_path.name),
        *([("Sheet", analysis.source.sheet)] if n_sheets > 1 else []),
        ("Wells analysed", f"{len(wells)} of {len(analysis.source.wells)}"),
        ("Analysis time range (h)", f"{analysis.recording.first_hour:g} \u2013 "
                                    f"{analysis.recording.last_hour:g}"),
        ("Time points", analysis.recording.n_points),
        ("Sampling interval (h)", analysis.time_interval),
        ("Outlier removal", f"running median {odd(s.outlier_window)} points, "
                            f"> {s.outlier_threshold:g} robust SD; "
                            f"{sum(analysis.outliers.values())} points removed"
         if s.remove_outliers else "Off"),
        ("Detrending", s.detrending_method),
        ("Trend parameter", trend.summary),
        ("Smoothing for peak detection", analysis.smooth_label),
        ("Peak min. separation (h)", s.peak_min_separation_hours),
        ("Peak min. prominence", s.peak_min_prominence or "Off"),
        ("Peaks and troughs", "Edited manually" if edited else "Automatic detection"),
        ("Actogram period (h)", periods[wells[0]] if len(set(periods.values())) == 1
         else "Set per well (see the regression sheet)"),
        ("Regression", "Linear regression of the selected peak/trough times on cycle number; "
                       "phase = fitted event time modulo the period, counted from 0 h"),
        ("Sine fit model", s.fit_model),
        ("Sine fit range (h)", f"{fit_start:g} \u2013 {fit_end:g}"),
        ("Sine fit period bounds (h)", f"{s.fit_min_period:g} \u2013 {s.fit_max_period:g}"),
        ("Sine fit definitions",
         "y = A\u00b7exp(-d\u00b7t')\u00b7sin(2\u03c0t'/P + \u03c6) + C with t' = time since the "
         "start of the fit range; A = amplitude at the start of the range; "
         "Peak phase = time of a fitted peak modulo P, counted from 0 h"),
        ("Processed (UTC)", utc_now_string()),
    ])

    report(0.27, "Writing the results workbook...")
    with pd.ExcelWriter(results_xlsx, engine="openpyxl") as writer:
        for name, frame in (("Note", note), ("Summary", summary), ("Sine Fit", fit_table),
                            ("Period-Phase Regression", reg_table),
                            ("Regression Points", reg_points), ("Peaks and Troughs", table)):
            frame.to_excel(writer, sheet_name=name, index=False)
            autosize(writer, name, frame)
    print(f"Written: {results_xlsx.name}")

    report(0.32, "Writing the data workbook...")
    with pd.ExcelWriter(data_xlsx, engine="openpyxl") as writer:
        note.to_excel(writer, sheet_name="Note", index=False)
        raw_out = raw.copy()
        raw_out.insert(1, "Temperature (\u00b0C)", analysis.recording.temperature.to_numpy())
        for name, frame in (("Raw Data", raw_out), ("Smoothed Raw", analysis.smoothed_raw),
                            (sheet_name(f"Trend ({trend.short_label})"), trend.data),
                            ("Detrended", detrended), ("Detrended Smoothed", smoothed)):
            for_excel(frame).to_excel(writer, sheet_name=name, index=False)
        if s.well_sheets:
            for i, well in enumerate(wells):
                per_well = pd.DataFrame({
                    TIME_OUT: raw[TIME], "Raw": raw[well],
                    "Smoothed raw": analysis.smoothed_raw[well], "Trend": trend.data[well],
                    "Detrended": detrended[well], "Detrended smoothed": smoothed[well],
                })
                detail = fit_details[well]
                if detail["fit"] is not None:
                    curve = pd.Series(np.nan, index=raw.index)
                    curve.loc[detail["times"].index] = fit_curve(detail["fit"], detail["times"])
                    per_well["Sine fit (detrended)"] = curve
                per_well.to_excel(writer, sheet_name=sheet_name(f"Well {well}"), index=False)
                if i % 10 == 0:
                    report(0.32 + 0.1 * i / len(wells), f"Writing sheet for well {well}...")
    print(f"Written: {data_xlsx.name}")

    save_peaks_json(peaks_json, peaks, smoothed[TIME], analysis.input_path.name,
                    {w: regressions[w]["used"] for w in wells}, periods)
    print(f"Written: {peaks_json.name}")

    report(0.45, "Plotting...")
    print("\nPlotting...")
    axis = make_axis_style(raw[TIME], MAJOR_TICK_HOURS[s.major_ticks], MINOR_TICK_HOURS[s.minor_ticks])
    rec = analysis.recording
    with PdfPages(plots_pdf) as pdf:
        plot_overview(pdf, title=f"{title} \u2013 raw data", wells=wells, layout=s.overview_layout,
                      scatter=raw, line=trend.data, line_label=f"trend ({s.detrending_method})",
                      y_label=s.value_label, axis=axis)
        plot_overview(pdf, title=f"{title} \u2013 detrended data", wells=wells,
                      layout=s.overview_layout, scatter=detrended, line=smoothed,
                      line_label=analysis.smooth_label, y_label="Detrended value", axis=axis)
        plot_temperature_page(pdf, rec, title, axis)
        plot_plate_map(pdf, analysis.source.wells, wells, title)
        plot_summary_pages(pdf, summary, title, s)
        if s.well_pages:
            for i, well in enumerate(wells):
                report(0.5 + 0.3 * i / len(wells), f"Plotting well {well}...")
                plot_well_page(
                    pdf, well=well, title=title, value_label=s.value_label, raw=raw, trend=trend,
                    detrended=detrended, smoothed=smoothed, axis=axis,
                    smooth_label=analysis.smooth_label, peaks_troughs=peaks[well],
                    reg_selection=regressions[well]["used"], regression=regressions[well],
                    day_length=periods[well], first_hour=rec.first_hour, last_hour=rec.last_hour,
                    label_peaks=s.label_detrended_plot, label_actogram=s.label_actogram,
                )
        if s.fit_pages:
            for i, well in enumerate(wells):
                report(0.8 + 0.18 * i / len(wells), f"Plotting fit for well {well}...")
                if fit_details[well]["fit"] is not None:
                    plot_fit_page(pdf, well=well, title=title, value_label=s.value_label, raw=raw,
                                  trend=trend, detrended=detrended, detail=fit_details[well],
                                  axis=axis)
    print(f"Written: {plots_pdf.name}")

    elapsed = (utc_now() - start_time).total_seconds()
    print(f"\nFinished in {elapsed:.0f} s at {utc_now_string()} (UTC)")
    report(1.0, "Done")
    return [results_xlsx, data_xlsx, plots_pdf, peaks_json]


def run_analysis(input_path, output_dir, settings: Settings, progress=None) -> list[Path]:
    """analyze() + write_outputs() without the review step (for scripts)."""
    analysis = analyze(input_path, settings, progress)
    return write_outputs(analysis, output_dir, progress=progress)


# =========================================================================== #
# 2. Desktop GUI
# =========================================================================== #

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except ImportError:  # the processing core above stays usable from scripts
    tk = None


def _fit_geometry(window, width: int, height: int) -> None:
    """Set the window size, but never larger than the screen."""
    width = min(width, window.winfo_screenwidth() - 20)
    height = min(height, window.winfo_screenheight() - 80)
    window.geometry(f"{width}x{height}")


class _QueueWriter:
    """File-like object that forwards print() output to the GUI thread."""

    def __init__(self, q: queue.Queue) -> None:
        self._q = q

    def write(self, text: str) -> int:
        if text:
            self._q.put(("log", text))
        return len(text)

    def flush(self) -> None:
        pass


class App:
    PREVIEW_MAX_LEGEND = 12

    def __init__(self, root) -> None:
        self.root = root
        self.queue: queue.Queue = queue.Queue()
        self.worker: threading.Thread | None = None
        self.source: Recording | None = None
        self.source_key: tuple | None = None
        self.analysis: Analysis | None = None
        self.pending_out_dir = ""
        self.last_output_dir: Path | None = None
        self._spin_specs: list = []
        self._name_from_file = True

        root.title(f"{APP_TITLE} {APP_VERSION}")
        _fit_geometry(root, 1320, 900)
        root.minsize(1100, 700)
        self._make_variables()
        self._build_ui()
        self._on_method_change()
        self._preview_job = None
        for var in (self.start_var, self.end_var, self.fit_start_var, self.fit_end_var,
                    self.value_label_var):
            var.trace_add("write", lambda *_: self._schedule_preview())
        self.root.after(100, self._poll_queue)

    def _schedule_preview(self) -> None:
        """Redraw the preview shortly after the last change (typing stays smooth)."""
        if self._preview_job is not None:
            self.root.after_cancel(self._preview_job)
        self._preview_job = self.root.after(400, self._run_preview_job)

    def _run_preview_job(self) -> None:
        self._preview_job = None
        self._draw_preview()

    # ------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------ #

    def _make_variables(self) -> None:
        d = Settings()
        self.in_var = tk.StringVar()
        self.sheet_var = tk.StringVar()
        self.out_var = tk.StringVar()
        self.summary_var = tk.StringVar(value="Select an Excel file (.xlsx) to begin.")
        self.name_var = tk.StringVar(value=d.experiment_name)
        self.title_var = tk.StringVar(value=d.experiment_title)
        self.value_label_var = tk.StringVar(value=d.value_label)
        self.pattern_var = tk.StringVar()
        self.wells_info_var = tk.StringVar(value="")

        self.start_var = tk.DoubleVar(value=d.start_hour)
        self.end_var = tk.DoubleVar(value=d.end_hour)
        self.outlier_var = tk.BooleanVar(value=d.remove_outliers)
        self.outlier_window_var = tk.IntVar(value=d.outlier_window)
        self.outlier_thr_var = tk.DoubleVar(value=d.outlier_threshold)

        self.method_var = tk.StringVar(value=d.detrending_method)
        self.cutoff_var = tk.DoubleVar(value=d.sinc_cutoff_hours)
        self.order_var = tk.IntVar(value=d.sinc_order)
        self.ma_var = tk.DoubleVar(value=d.ma_window_hours)
        self.poly_var = tk.IntVar(value=d.poly_degree)

        self.smooth_var = tk.IntVar(value=d.smoothing_points)
        self.sep_var = tk.DoubleVar(value=d.peak_min_separation_hours)
        self.prom_var = tk.DoubleVar(value=d.peak_min_prominence)

        self.layout_var = tk.StringVar(value=d.overview_layout)
        self.major_var = tk.StringVar(value=d.major_ticks)
        self.minor_var = tk.StringVar(value=d.minor_ticks)
        self.label_det_var = tk.BooleanVar(value=d.label_detrended_plot)
        self.label_act_var = tk.BooleanVar(value=d.label_actogram)
        self.acto_var = tk.DoubleVar(value=d.actogram_period)
        self.well_pages_var = tk.BooleanVar(value=d.well_pages)
        self.fit_pages_var = tk.BooleanVar(value=d.fit_pages)
        self.well_sheets_var = tk.BooleanVar(value=d.well_sheets)

        self.model_var = tk.StringVar(value=d.fit_model)
        self.fit_start_var = tk.DoubleVar(value=d.fit_start_hour)
        self.fit_end_var = tk.DoubleVar(value=d.fit_end_hour)
        self.pmin_var = tk.DoubleVar(value=d.fit_min_period)
        self.pmax_var = tk.DoubleVar(value=d.fit_max_period)

        self.review_var = tk.BooleanVar(value=True)
        self.status_var = tk.StringVar(value="Ready.")

    # ------------------------------------------------------------------ #
    # UI construction helpers
    # ------------------------------------------------------------------ #

    @staticmethod
    def _box(parent, text: str) -> ttk.LabelFrame:
        box = ttk.LabelFrame(parent, text=text, padding=8)
        box.pack(fill="x", pady=(0, 8))
        box.columnconfigure(1, weight=1)
        return box

    @staticmethod
    def _label(box, row: int, text: str) -> ttk.Label:
        label = ttk.Label(box, text=text)
        label.grid(row=row, column=0, sticky="w", padx=(0, 10), pady=2)
        return label

    def _entry(self, box, row: int, text: str, var) -> ttk.Entry:
        self._label(box, row, text)
        entry = ttk.Entry(box, textvariable=var)
        entry.grid(row=row, column=1, columnspan=2, sticky="ew", pady=2)
        return entry

    def _combo(self, box, row: int, text: str, var, values) -> ttk.Combobox:
        self._label(box, row, text)
        combo = ttk.Combobox(box, textvariable=var, values=list(values), state="readonly", width=30)
        combo.grid(row=row, column=1, columnspan=2, sticky="ew", pady=2)
        return combo

    def _spin(self, box, row: int, text: str, var, lo, hi, step, hint: str = ""):
        label = self._label(box, row, text)
        is_float = isinstance(step, float)
        spin = ttk.Spinbox(box, from_=lo, to=hi, increment=step, textvariable=var, width=9,
                           justify="right", format="%.2f" if is_float else "%.0f")
        spin.grid(row=row, column=1, sticky="w", pady=2)
        if hint:
            ttk.Label(box, text=hint, foreground="#666666").grid(row=row, column=2, sticky="w",
                                                                 padx=(6, 0))
        self._spin_specs.append((text, var, lo, hi, spin))
        return label, spin

    def _check(self, box, row: int, text: str, var, command=None) -> ttk.Checkbutton:
        check = ttk.Checkbutton(box, text=text, variable=var, command=command)
        check.grid(row=row, column=0, columnspan=3, sticky="w", pady=1)
        return check

    def _validate_spinners(self) -> None:
        for text, var, lo, hi, spin in self._spin_specs:
            if spin.instate(["disabled"]):
                continue
            try:
                value = float(var.get())
            except (tk.TclError, ValueError):
                raise ValueError(f"'{text}' must be a number.") from None
            if not lo <= value <= hi:
                raise ValueError(f"'{text}' must be between {lo} and {hi}.")

    # ------------------------------------------------------------------ #
    # UI construction
    # ------------------------------------------------------------------ #

    def _build_ui(self) -> None:
        root = self.root
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        left = ttk.Frame(root, padding=(10, 10, 4, 10))
        left.grid(row=0, column=0, sticky="nsew")
        right = ttk.Frame(root, padding=(4, 10, 10, 10))
        right.grid(row=0, column=1, sticky="nsew")

        tabs = ttk.Notebook(left, width=440)
        tabs.pack(fill="both", expand=True)
        data_tab = ttk.Frame(tabs, padding=8)
        proc_tab = ttk.Frame(tabs, padding=8)
        out_tab = ttk.Frame(tabs, padding=8)
        tabs.add(data_tab, text=" 1. Data ")
        tabs.add(proc_tab, text=" 2. Processing ")
        tabs.add(out_tab, text=" 3. Fit && output ")

        self._build_files(data_tab)
        self._build_wells(data_tab)
        self._build_range(proc_tab)
        self._build_detrending(proc_tab)
        self._build_peaks(proc_tab)
        self._build_fit(out_tab)
        self._build_plots(out_tab)
        self._build_right(right)

    def _build_files(self, parent) -> None:
        box = self._box(parent, "Files")
        self._label(box, 0, "Excel file")
        ttk.Entry(box, textvariable=self.in_var, width=30).grid(row=0, column=1, sticky="ew", pady=2)
        ttk.Button(box, text="Browse...", command=self._browse_input).grid(row=0, column=2, padx=(6, 0))
        self._label(box, 1, "Sheet")
        self.sheet_combo = ttk.Combobox(box, textvariable=self.sheet_var, state="readonly", width=28)
        self.sheet_combo.grid(row=1, column=1, sticky="ew", pady=2)
        self.sheet_combo.bind("<<ComboboxSelected>>", lambda _e: self._load_source())
        ttk.Button(box, text="Reload", command=self._load_source).grid(row=1, column=2, padx=(6, 0))
        self._label(box, 2, "Output folder")
        ttk.Entry(box, textvariable=self.out_var, width=30).grid(row=2, column=1, sticky="ew", pady=2)
        ttk.Button(box, text="Browse...", command=self._browse_output).grid(row=2, column=2, padx=(6, 0))
        ttk.Label(box, textvariable=self.summary_var, wraplength=400, foreground="#444444").grid(
            row=3, column=0, columnspan=3, sticky="w", pady=(4, 0))

        box = self._box(parent, "Experiment")
        name_entry = self._entry(box, 0, "Name (file prefix)", self.name_var)
        name_entry.bind("<Key>", lambda _e: setattr(self, "_name_from_file", False))
        self._entry(box, 1, "Title (optional)", self.title_var)
        self._entry(box, 2, "Y-axis label", self.value_label_var)

    def _build_wells(self, parent) -> None:
        box = ttk.LabelFrame(parent, text="Wells to analyse (selected = analysed)", padding=8)
        box.pack(fill="both", expand=True)
        box.columnconfigure(0, weight=1)
        box.rowconfigure(0, weight=1)
        frame = ttk.Frame(box)
        frame.grid(row=0, column=0, columnspan=4, sticky="nsew")
        self.well_list = tk.Listbox(frame, selectmode="extended", exportselection=False, height=10)
        scroll = ttk.Scrollbar(frame, orient="vertical", command=self.well_list.yview)
        self.well_list.configure(yscrollcommand=scroll.set)
        self.well_list.pack(side="left", fill="both", expand=True)
        scroll.pack(side="left", fill="y")
        self.well_list.bind("<<ListboxSelect>>", lambda _e: self._on_wells_changed())

        ttk.Button(box, text="All", width=6, command=lambda: self._select_wells(True)).grid(
            row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Button(box, text="None", width=6, command=lambda: self._select_wells(False)).grid(
            row=1, column=1, sticky="w", pady=(6, 0))
        ttk.Label(box, textvariable=self.wells_info_var, foreground="#444444").grid(
            row=1, column=2, columnspan=2, sticky="e", pady=(6, 0))
        pattern = ttk.Frame(box)
        pattern.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(6, 0))
        ttk.Label(pattern, text="Pattern").pack(side="left")
        entry = ttk.Entry(pattern, textvariable=self.pattern_var, width=16)
        entry.pack(side="left", padx=4, fill="x", expand=True)
        entry.bind("<Return>", lambda _e: self._select_pattern(replace_selection=True))
        ttk.Button(pattern, text="Select", command=lambda: self._select_pattern(True)).pack(side="left")
        ttk.Button(pattern, text="Add", command=lambda: self._select_pattern(False)).pack(
            side="left", padx=(4, 0))
        ttk.Label(box, text="e.g.  K*   D4 D5   row:J   col:5   (comma/space separated)",
                  foreground="#666666").grid(row=3, column=0, columnspan=4, sticky="w")

    def _build_range(self, parent) -> None:
        box = self._box(parent, "Time range and outliers")
        self._spin(box, 0, "Analyse from (h)", self.start_var, 0, 100000, 0.25)
        self._spin(box, 1, "Analyse to (h)", self.end_var, 0, 100000, 0.25, "0 = to the end")
        self._check(box, 2, "Remove spikes (running-median outlier filter)", self.outlier_var,
                    self._on_outlier_toggle)
        self._outlier_widgets = [
            self._spin(box, 3, "   Window (points)", self.outlier_window_var, 3, 101, 2),
            self._spin(box, 4, "   Threshold (robust SD)", self.outlier_thr_var, 2, 50, 0.5),
        ]
        ttk.Label(box, text="Tip: exclude an initial transient with 'Analyse from'.",
                  foreground="#666666").grid(row=5, column=0, columnspan=3, sticky="w")
        self._on_outlier_toggle()

    def _build_detrending(self, parent) -> None:
        box = self._box(parent, "Detrending")
        combo = self._combo(box, 0, "Method", self.method_var, DETRENDING_METHODS)
        combo.bind("<<ComboboxSelected>>", lambda _e: self._on_method_change())
        self._sinc_widgets = [
            self._spin(box, 1, "Sinc cutoff period (h)", self.cutoff_var, 4, 1000, 1.0),
            self._spin(box, 2, "Sinc filter order", self.order_var, 0, 100001, 2, "0 = automatic"),
        ]
        self._ma_widgets = [self._spin(box, 3, "Moving-average window (h)", self.ma_var, 1, 1000, 1.0)]
        self._poly_widgets = [self._spin(box, 4, "Polynomial degree", self.poly_var, 1, 8, 1)]

    def _build_peaks(self, parent) -> None:
        box = self._box(parent, "Peak / trough detection")
        self._spin(box, 0, "Smoothing (points, odd)", self.smooth_var, 1, 301, 2)
        self._spin(box, 1, "Min. separation (h)", self.sep_var, 0, 200, 0.5)
        self._spin(box, 2, "Min. prominence", self.prom_var, 0, 1e9, 0.1, "0 = off")
        ttk.Label(box, text="Detection runs on the smoothed detrended data.",
                  foreground="#666666").grid(row=3, column=0, columnspan=3, sticky="w")

    def _build_fit(self, parent) -> None:
        box = self._box(parent, "Sine-curve fitting")
        self._combo(box, 0, "Model", self.model_var, FIT_MODELS)
        self._spin(box, 1, "Fit from (h)", self.fit_start_var, 0, 100000, 0.25)
        self._spin(box, 2, "Fit to (h)", self.fit_end_var, 0, 100000, 0.25, "0 = to the end")
        self._spin(box, 3, "Min. period (h)", self.pmin_var, 1, 1000, 0.5)
        self._spin(box, 4, "Max. period (h)", self.pmax_var, 1, 1000, 0.5)

    def _build_plots(self, parent) -> None:
        box = self._box(parent, "Plots and output")
        self._combo(box, 0, "Overview layout", self.layout_var, OVERVIEW_LAYOUTS)
        self._combo(box, 1, "Major ticks", self.major_var, MAJOR_TICK_HOURS)
        self._combo(box, 2, "Minor ticks", self.minor_var, MINOR_TICK_HOURS)
        self._spin(box, 3, "Actogram period (h)", self.acto_var, ACTOGRAM_MIN_PERIOD_HOURS,
                   ACTOGRAM_MAX_PERIOD_HOURS, 0.1, "default for all wells")
        self._check(box, 4, "Label peak/trough times in detrended plots", self.label_det_var)
        self._check(box, 5, "Label peak/trough times in actograms", self.label_act_var)
        self._check(box, 6, "PDF: one page per well (data + actogram)", self.well_pages_var)
        self._check(box, 7, "PDF: one sine-fit page per well", self.fit_pages_var)
        self._check(box, 8, "Data workbook: one sheet per well", self.well_sheets_var)

    def _build_right(self, parent) -> None:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure

        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(3, weight=1)
        buttons = ttk.Frame(parent)
        buttons.grid(row=0, column=0, sticky="ew")
        self.run_button = ttk.Button(buttons, text="Analyze \u25b6", command=self._on_run)
        self.run_button.pack(side="left")
        ttk.Checkbutton(buttons, text="Review peaks/troughs before export",
                        variable=self.review_var).pack(side="left", padx=(10, 0))
        self.open_button = ttk.Button(buttons, text="Open output folder",
                                      command=self._open_output, state="disabled")
        self.open_button.pack(side="right")

        self.progress = ttk.Progressbar(parent, mode="determinate", maximum=100)
        self.progress.grid(row=1, column=0, sticky="ew", pady=(10, 2))
        ttk.Label(parent, textvariable=self.status_var).grid(row=2, column=0, sticky="w")

        self.right_tabs = ttk.Notebook(parent)
        self.right_tabs.grid(row=3, column=0, sticky="nsew", pady=(8, 0))
        preview = ttk.Frame(self.right_tabs)
        log_frame = ttk.Frame(self.right_tabs)
        self.right_tabs.add(preview, text=" Preview ")
        self.right_tabs.add(log_frame, text=" Log ")

        bar = ttk.Frame(preview)
        bar.pack(fill="x", pady=(4, 2))
        ttk.Button(bar, text="Refresh preview", command=self._draw_preview).pack(side="left")
        ttk.Label(bar, text="  Shows the selected wells; grey band = outside the analysis range, "
                            "green lines = sine-fit range.", foreground="#555555").pack(side="left")
        self.preview_fig = Figure(figsize=(8, 6), dpi=100)
        grid = self.preview_fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
        self.pv_ax = self.preview_fig.add_subplot(grid[0])
        self.pv_temp = self.preview_fig.add_subplot(grid[1], sharex=self.pv_ax)
        self.preview_canvas = FigureCanvasTkAgg(self.preview_fig, master=preview)
        toolbar = NavigationToolbar2Tk(self.preview_canvas, preview, pack_toolbar=False)
        toolbar.pack(side="bottom", fill="x")
        self.preview_canvas.get_tk_widget().pack(side="top", fill="both", expand=True)

        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        self.log = tk.Text(log_frame, wrap="word", state="disabled", height=20)
        scroll = ttk.Scrollbar(log_frame, orient="vertical", command=self.log.yview)
        self.log.configure(yscrollcommand=scroll.set)
        self.log.grid(row=0, column=0, sticky="nsew")
        scroll.grid(row=0, column=1, sticky="ns")
        self._draw_preview()

    # ------------------------------------------------------------------ #
    # Settings-panel events
    # ------------------------------------------------------------------ #

    @staticmethod
    def _enable(group, enabled: bool) -> None:
        for label, spin in group:
            label.state(["!disabled"] if enabled else ["disabled"])
            spin.state(["!disabled"] if enabled else ["disabled"])

    def _on_method_change(self) -> None:
        method = self.method_var.get()
        self._enable(self._sinc_widgets, method == "Sinc Filter")
        self._enable(self._ma_widgets, method == "Moving Average")
        self._enable(self._poly_widgets, method == "Polynomial")

    def _on_outlier_toggle(self) -> None:
        self._enable(self._outlier_widgets, bool(self.outlier_var.get()))

    def _browse_input(self) -> None:
        path = filedialog.askopenfilename(
            title="Select the time-series Excel file",
            filetypes=[("Excel files", "*.xlsx *.xlsm *.xls"), ("All files", "*.*")],
        )
        if not path:
            return
        self.in_var.set(path)
        self.out_var.set(str(Path(path).parent))
        try:
            sheets = list_sheets(path)
        except Exception as error:  # noqa: BLE001 - reported to the user
            messagebox.showerror(APP_TITLE, f"Could not open the workbook:\n{error}")
            return
        self.sheet_combo.configure(values=sheets)
        self.sheet_var.set(sheets[0] if sheets else "")
        if self._name_from_file or not self.name_var.get().strip():
            self.name_var.set(safe_filename(Path(path).stem))
            self._name_from_file = True
        self._load_source()

    def _browse_output(self) -> None:
        path = filedialog.askdirectory(title="Select the output folder")
        if path:
            self.out_var.set(path)

    def _load_source(self) -> None:
        path = self.in_var.get().strip()
        if not path or not Path(path).is_file():
            return
        self.root.configure(cursor="watch")
        self.root.update_idletasks()
        try:
            source = read_timeseries_excel(path, self.sheet_var.get() or None)
        except Exception as error:  # noqa: BLE001
            self.source, self.source_key = None, None
            self.summary_var.set(f"Could not read this sheet: {error}")
            self.well_list.delete(0, "end")
            self._draw_preview()
            return
        finally:
            self.root.configure(cursor="")
        self.source = source
        self.source_key = (path, self.sheet_var.get())
        self.summary_var.set(source.describe())
        self.well_list.delete(0, "end")
        for well in source.wells:
            self.well_list.insert("end", well)
        self.well_list.selection_set(0, "end")
        self._on_wells_changed()

    def _selected_wells(self) -> list[str]:
        if self.source is None:
            return []
        wells = self.source.wells
        return [wells[i] for i in self.well_list.curselection()]

    def _select_wells(self, select: bool) -> None:
        if select:
            self.well_list.selection_set(0, "end")
        else:
            self.well_list.selection_clear(0, "end")
        self._on_wells_changed()

    def _select_pattern(self, replace_selection: bool) -> None:
        if self.source is None:
            return
        chosen = set(match_wells(self.source.wells, self.pattern_var.get()))
        if not chosen:
            self.status_var.set("No well matches that pattern.")
            return
        if replace_selection:
            self.well_list.selection_clear(0, "end")
        for i, well in enumerate(self.source.wells):
            if well in chosen:
                self.well_list.selection_set(i)
        self._on_wells_changed()

    def _on_wells_changed(self) -> None:
        total = len(self.source.wells) if self.source else 0
        self.wells_info_var.set(f"{len(self._selected_wells())} of {total} selected")
        self._draw_preview()

    def _draw_preview(self) -> None:
        ax, ax_t = self.pv_ax, self.pv_temp
        ax.clear()
        ax_t.clear()
        if self.source is None:
            ax.text(0.5, 0.5, "No data loaded", ha="center", va="center", transform=ax.transAxes,
                    color="gray")
            self.preview_canvas.draw_idle()
            return
        data = self.source.data
        wells = self._selected_wells()
        many = len(wells) > self.PREVIEW_MAX_LEGEND
        for well in wells:
            ax.plot(data[TIME], data[well], linewidth=0.5 if many else 0.8, alpha=0.6 if many else 0.9,
                    label=None if many else well)
        if wells and not many:
            ax.legend(fontsize=7, ncol=2, loc="upper right")
        ax.set_ylabel(self.value_label_var.get() or "Value")
        ax.set_title(f"{Path(self.source.filename).name}  \u2013 "
                     f"{len(wells)} well(s)", fontsize=9)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax_t.plot(data[TIME], self.source.temperature, color="tab:orange", linewidth=0.8)
        ax_t.set_ylabel("Temp. (\u00b0C)")
        ax_t.set_xlabel("Time (h)")
        ax_t.grid(True, linestyle="--", alpha=0.5)
        first, last = self.source.first_hour, self.source.last_hour
        try:
            start = float(self.start_var.get())
            end = float(self.end_var.get()) or last
            fit_start = float(self.fit_start_var.get())
            fit_end = float(self.fit_end_var.get()) or last
        except (tk.TclError, ValueError):
            start, end, fit_start, fit_end = first, last, first, last
        for a in (ax, ax_t):
            if start > first:
                a.axvspan(first, start, color="gray", alpha=0.2, lw=0)
            if end < last:
                a.axvspan(end, last, color="gray", alpha=0.2, lw=0)
        for x in (fit_start, fit_end):
            ax.axvline(x, color="green", linestyle=":", linewidth=1)
        ax.set_xlim(first, last)
        ax.xaxis.set_major_locator(MultipleLocator(24 if last - first > 48 else 6))
        ax.tick_params(labelbottom=False)
        self.preview_fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.09)
        self.preview_canvas.draw_idle()

    def _collect_settings(self) -> Settings:
        self._validate_spinners()
        name = self.name_var.get().strip()
        if not safe_filename(name):
            raise ValueError("Please enter an experiment name (used as the file prefix).")
        wells = self._selected_wells()
        if not wells:
            raise ValueError("Select at least one well.")
        start, end = float(self.start_var.get()), float(self.end_var.get())
        if end and end <= start + 24:
            raise ValueError("'Analyse to' must be 0 (the end) or at least 24 h after 'Analyse from'.")
        fit_start, fit_end = float(self.fit_start_var.get()), float(self.fit_end_var.get())
        if fit_end and fit_end <= fit_start + 12:
            raise ValueError("'Fit to' must be 0 (the end) or at least 12 h after 'Fit from'.")
        pmin, pmax = float(self.pmin_var.get()), float(self.pmax_var.get())
        if pmin >= pmax:
            raise ValueError("The minimum fit period must be smaller than the maximum.")
        return Settings(
            experiment_name=name, experiment_title=self.title_var.get().strip(),
            value_label=self.value_label_var.get().strip() or DEFAULT_VALUE_LABEL,
            sheet=self.sheet_var.get(), wells=wells, start_hour=start, end_hour=end,
            remove_outliers=bool(self.outlier_var.get()),
            outlier_window=int(self.outlier_window_var.get()),
            outlier_threshold=float(self.outlier_thr_var.get()),
            detrending_method=self.method_var.get(),
            sinc_cutoff_hours=float(self.cutoff_var.get()), sinc_order=int(self.order_var.get()),
            ma_window_hours=float(self.ma_var.get()), poly_degree=int(self.poly_var.get()),
            smoothing_points=odd(int(self.smooth_var.get())),
            peak_min_separation_hours=float(self.sep_var.get()),
            peak_min_prominence=float(self.prom_var.get()),
            overview_layout=self.layout_var.get(), major_ticks=self.major_var.get(),
            minor_ticks=self.minor_var.get(), label_detrended_plot=bool(self.label_det_var.get()),
            label_actogram=bool(self.label_act_var.get()),
            actogram_period=round(float(self.acto_var.get()), 2),
            well_pages=bool(self.well_pages_var.get()), fit_pages=bool(self.fit_pages_var.get()),
            well_sheets=bool(self.well_sheets_var.get()), fit_model=self.model_var.get(),
            fit_start_hour=fit_start, fit_end_hour=fit_end, fit_min_period=pmin, fit_max_period=pmax,
        )

    # ------------------------------------------------------------------ #
    # Running
    # ------------------------------------------------------------------ #

    def _on_run(self) -> None:
        if self.worker is not None and self.worker.is_alive():
            return
        in_path = self.in_var.get().strip()
        if not in_path or not Path(in_path).is_file():
            messagebox.showerror(APP_TITLE, "Please select an existing Excel file.")
            return
        if self.source is None or self.source_key != (in_path, self.sheet_var.get()):
            self._load_source()
            if self.source is None:
                messagebox.showerror(APP_TITLE, "The selected sheet could not be read.")
                return
        out_dir = self.out_var.get().strip() or str(Path(in_path).parent)
        self.out_var.set(out_dir)
        try:
            settings = self._collect_settings()
        except ValueError as error:
            messagebox.showerror(APP_TITLE, str(error))
            return
        existing = [p.name for p in expected_outputs(out_dir, settings) if p.exists()]
        if existing and not messagebox.askyesno(
            APP_TITLE,
            "These files already exist and will be overwritten:\n\n" + "\n".join(existing)
            + "\n\nContinue?",
        ):
            return
        self._clear_log()
        self.right_tabs.select(1)
        self.progress["value"] = 0
        self.status_var.set("Starting...")
        self.run_button.configure(state="disabled")
        self.open_button.configure(state="disabled")
        self.pending_out_dir = out_dir
        self.analysis = None
        self._run_in_thread(self._analysis_worker, in_path, settings, self.source)

    def _open_output(self) -> None:
        folder = self.last_output_dir
        if folder is None:
            return
        try:
            if sys.platform.startswith("win"):
                os.startfile(str(folder))  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.run(["open", str(folder)], check=False)
            else:
                subprocess.run(["xdg-open", str(folder)], check=False)
        except OSError as error:
            messagebox.showerror(APP_TITLE, f"Could not open the folder: {error}")

    def _run_in_thread(self, target, *args) -> None:
        self.worker = threading.Thread(target=target, args=args, daemon=True)
        self.worker.start()

    def _progress_callback(self):
        return lambda fraction, text: self.queue.put(("progress", fraction, text))

    def _analysis_worker(self, in_path, settings, source) -> None:
        writer = _QueueWriter(self.queue)
        try:
            with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
                analysis = analyze(in_path, settings, self._progress_callback(), source=source)
            self.queue.put(("analyzed", analysis))
        except Exception as error:  # noqa: BLE001
            self.queue.put(("error", traceback.format_exc(), str(error)))

    def _export_worker(self, analysis, out_dir, peaks, reg_selection, periods) -> None:
        writer = _QueueWriter(self.queue)
        try:
            with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
                outputs = write_outputs(analysis, out_dir, peaks, self._progress_callback(),
                                        reg_selection, periods)
            self.queue.put(("done", outputs))
        except Exception as error:  # noqa: BLE001
            self.queue.put(("error", traceback.format_exc(), str(error)))

    def _on_analyzed(self, analysis: Analysis) -> None:
        self.analysis = analysis
        if self.review_var.get():
            self.status_var.set("Review the peaks and troughs...")
            ReviewWindow(self.root, analysis, on_ok=self._start_export,
                         on_cancel=self._on_review_cancel)
        else:
            self._start_export(analysis.auto_peaks, default_reg_selection(analysis.auto_peaks),
                               analysis.settings.actogram_period)

    def _start_export(self, peaks, reg_selection, periods) -> None:
        self.progress["value"] = 0
        self.status_var.set("Writing output files...")
        self._run_in_thread(self._export_worker, self.analysis, self.pending_out_dir,
                            peaks, reg_selection, periods)

    def _on_review_cancel(self) -> None:
        self.run_button.configure(state="normal")
        self.progress["value"] = 0
        self.status_var.set("Cancelled. No files were written.")
        self._append_log("\nCancelled in the review window. No files were written.\n")

    def _poll_queue(self) -> None:
        try:
            while True:
                message = self.queue.get_nowait()
                kind = message[0]
                if kind == "log":
                    self._append_log(message[1])
                elif kind == "progress":
                    self.progress["value"] = message[1] * 100
                    self.status_var.set(message[2])
                elif kind == "analyzed":
                    self._on_analyzed(message[1])
                elif kind == "done":
                    self._on_done(message[1])
                elif kind == "error":
                    self._on_error(message[1], message[2])
        except queue.Empty:
            pass
        self.root.after(100, self._poll_queue)

    def _on_done(self, outputs) -> None:
        self.run_button.configure(state="normal")
        self.progress["value"] = 100
        self.status_var.set("Done.")
        if outputs:
            self.last_output_dir = Path(outputs[0]).parent
            self.open_button.configure(state="normal")
        names = "\n".join(Path(p).name for p in outputs)
        messagebox.showinfo(APP_TITLE, f"Finished. Files written to:\n{self.last_output_dir}\n\n{names}")

    def _on_error(self, details: str, short: str) -> None:
        self.run_button.configure(state="normal")
        self.status_var.set("Failed.")
        self._append_log("\n" + details)
        messagebox.showerror(APP_TITLE, f"The analysis failed:\n\n{short}")

    def _append_log(self, text: str) -> None:
        self.log.configure(state="normal")
        self.log.insert("end", text)
        self.log.see("end")
        self.log.configure(state="disabled")

    def _clear_log(self) -> None:
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")


class ReviewWindow:
    """Modal window for checking and hand-editing the peaks and troughs.

    Top plot     : raw data + trend (context).
    Middle plot  : detrended data; add / remove peaks and troughs; optional sine fit.
    Bottom plot  : double-plotted actogram with an adjustable period; the points
                   picked with the mouse are used for the period/phase regression.
    """

    HIT_RADIUS_PX = 25
    ZOOM_WINDOW_HOURS = 24.0     # width of the view when zooming onto one peak/trough
    ZOOM_MIN_WIDTH_HOURS = 2.0   # the view is never narrower than this
    WHEEL_ZOOM_FACTOR = 0.85     # view width is multiplied by this per wheel step (in)
    WHEEL_PAN_FRACTION = 0.15    # Shift+wheel scrolls by this fraction of the view
    DRAG_THRESHOLD_PX = 6

    HELP_TEXT = (
        "Detrended plot (middle)\n"
        "  \u2022 'Remove': click a marker to delete that peak/trough.\n"
        "  \u2022 'Add peak' / 'Add trough': click near the curve; the point snaps to the local "
        "extreme of the smoothed curve within \u00b1 1 h.\n"
        "  \u2022 'Show sine fit' overlays the fitted curve.\n\n"
        "Zooming the detrended plot (middle)\n"
        "  \u2022 Mouse wheel: zoom in/out around the pointer.  Shift+wheel: scroll in time.\n"
        "  \u2022 Right-click: zoom onto the nearest peak/trough.  Right-drag: pan.\n"
        "  \u2022 Middle-click or right double-click: full view (also 'Full view', Esc, "
        "and the toolbar's Home button).\n"
        "  \u2022 'Prev/Next marker' (or the , and . keys) step through the peaks and troughs "
        "one by one; the current one is ringed in orange.\n"
        "  \u2022 The raw-data plot (top) is an overview: the shaded band shows the zoomed "
        "range, and a click there moves the view to that time.\n"
        "  \u2022 Adding/removing markers works at any zoom level.\n\n"
        "Actogram (bottom)\n"
        "  \u2022 Click a point to use / not use it in the period-phase regression "
        "(ringed = used).\n"
        "  \u2022 Drag a box to select the points inside; right-drag to deselect them.\n"
        "  \u2022 The actogram period is set per well ('24 h', '= fit period', "
        "'Apply to all wells').\n\n"
        "General\n"
        "  \u2022 While a zoom/pan tool of the toolbar is active, clicks do not edit.\n"
        "  \u2022 Switch wells with the list or the Left/Right arrow keys.\n"
        "  \u2022 'Save edits...' / 'Load edits...' store the selection as a JSON file."
    )

    def __init__(self, parent, analysis: Analysis, on_ok, on_cancel) -> None:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure

        self.analysis = analysis
        self.on_ok = on_ok
        self.on_cancel = on_cancel
        self.wells = list(analysis.wells)
        self.data = analysis.detrended
        self.smooth = analysis.smoothed
        self.auto = analysis.auto_peaks
        self.selection = {w: (list(p), list(t)) for w, (p, t) in self.auto.items()}
        self.reg_sel = self._all_selected()
        default_period = self._clamp_period(float(analysis.settings.actogram_period))
        self.periods = {w: default_period for w in self.wells}
        self.well = self.wells[0]
        self._fit_cache: dict = {}
        self._marker_artists: list = []
        self._drag: dict | None = None
        self._drag_patch = None
        self._pan: dict | None = None        # right-button drag in the middle plot
        self._focus_row: int | None = None   # peak/trough the zoomed view is centred on
        self._view_span = None               # shaded band in the raw-data plot
        self._full_xlim = (0.0, 1.0)
        self._full_ylim = (0.0, 1.0)
        self._closed = False

        win = self.win = tk.Toplevel(parent)
        win.title("Review peaks and troughs")
        _fit_geometry(win, 1280, 1040)
        win.minsize(950, 680)
        # Not win.transient(parent): a transient (dialog) window gets only a Close button
        # on Windows and many Linux desktops.  As a normal window it also has Minimize and
        # Maximize; _link_to_parent() keeps it together with the (blocked) main window.
        win.resizable(True, True)
        win.protocol("WM_DELETE_WINDOW", self._cancel)
        self.parent = parent
        self._syncing_state = False
        win.columnconfigure(1, weight=1)
        win.rowconfigure(1, weight=1)

        self.mode_var = tk.StringVar(value="remove")
        self.status_var = tk.StringVar(value="")
        self.period_var = tk.StringVar(value=f"{self.day_length:g}")
        self.reg_var = tk.StringVar(value="")
        self.fit_var = tk.BooleanVar(value=True)

        # (the usage notes are behind the "Help" button, so the plots get the height)
        ttk.Frame(win, height=6).grid(row=0, column=0, columnspan=2)

        left = ttk.Frame(win, padding=(10, 0, 4, 0))
        left.grid(row=1, column=0, sticky="ns")
        ttk.Label(left, text="Well  (P peaks, T troughs, * edited)").pack(anchor="w")
        list_frame = ttk.Frame(left)
        list_frame.pack(fill="y", expand=True)
        self.listbox = tk.Listbox(list_frame, width=24, height=30, exportselection=False)
        scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.listbox.yview)
        self.listbox.configure(yscrollcommand=scrollbar.set)
        self.listbox.pack(side="left", fill="y")
        scrollbar.pack(side="left", fill="y")
        for well in self.wells:
            self.listbox.insert("end", self._list_text(well))
        self.listbox.selection_set(0)
        self.listbox.bind("<<ListboxSelect>>", self._on_list_select)

        centre = ttk.Frame(win, padding=(4, 0, 10, 0))
        centre.grid(row=1, column=1, sticky="nsew")

        controls = ttk.Frame(centre)
        controls.pack(fill="x", pady=(0, 4))
        ttk.Label(controls, text="Middle plot:").pack(side="left", padx=(0, 6))
        for text, value in (("Remove", "remove"), ("Add peak", "add_peak"),
                            ("Add trough", "add_trough")):
            ttk.Radiobutton(controls, text=text, value=value, variable=self.mode_var).pack(
                side="left", padx=(0, 12))
        ttk.Checkbutton(controls, text="Show sine fit", variable=self.fit_var,
                        command=lambda: self._draw_full(keep_view=True)).pack(
            side="left", padx=(6, 0))
        ttk.Button(controls, text="Help", width=6, command=self._show_help).pack(side="right")
        ttk.Button(controls, text="Next \u25b6", command=lambda: self._step(1)).pack(
            side="right", padx=(0, 12))
        ttk.Button(controls, text="\u25c0 Previous", command=lambda: self._step(-1)).pack(
            side="right", padx=4)
        ttk.Button(controls, text="Reset this well", command=self._reset_well).pack(
            side="right", padx=(0, 12))

        zoom = ttk.Frame(centre)
        zoom.pack(fill="x", pady=(0, 4))
        ttk.Label(zoom, text="Zoom:").pack(side="left", padx=(0, 6))
        ttk.Button(zoom, text="\u25c0 Prev marker", command=lambda: self._jump_marker(-1)).pack(
            side="left")
        ttk.Button(zoom, text="Next marker \u25b6", command=lambda: self._jump_marker(1)).pack(
            side="left", padx=(4, 0))
        ttk.Button(zoom, text="Full view", command=self._full_view).pack(side="left", padx=(4, 0))
        ttk.Label(zoom, text="wheel: zoom   Shift+wheel: scroll   right-click: zoom to marker   "
                             "right-drag: pan   middle-click: full view",
                  foreground="gray40").pack(side="left", padx=(12, 0))

        acto = ttk.Frame(centre)
        acto.pack(fill="x", pady=(0, 4))
        ttk.Label(acto, text="Actogram period (h), this well:").pack(side="left", padx=(0, 4))
        self.period_spin = ttk.Spinbox(
            acto, from_=ACTOGRAM_MIN_PERIOD_HOURS, to=ACTOGRAM_MAX_PERIOD_HOURS, increment=0.1,
            textvariable=self.period_var, width=7, justify="right", format="%.1f",
            command=self._apply_period)
        self.period_spin.pack(side="left")
        self.period_spin.bind("<Return>", self._apply_period)
        self.period_spin.bind("<FocusOut>", self._apply_period)
        ttk.Button(acto, text="24 h", width=5, command=lambda: self._set_period(24.0)).pack(
            side="left", padx=(6, 0))
        ttk.Button(acto, text="= fit period", command=self._period_from_fit).pack(
            side="left", padx=(6, 0))
        ttk.Button(acto, text="Apply to all wells", command=self._period_to_all).pack(
            side="left", padx=(6, 0))
        ttk.Separator(acto, orient="vertical").pack(side="left", fill="y", padx=12)
        ttk.Label(acto, text="Regression points:").pack(side="left", padx=(0, 6))
        ttk.Button(acto, text="Select all", command=lambda: self._select_all(True)).pack(side="left")
        ttk.Button(acto, text="Clear", command=lambda: self._select_all(False)).pack(
            side="left", padx=(4, 0))

        # Constrained layout sizes each gap to what is in it: the raw-data graph has no
        # x tick labels, so it sits close above the detrended graph, while the actogram
        # still gets room for the axis label above it.  Outside legends are included.
        self.fig = Figure(figsize=(8, 8), dpi=100, layout="constrained")
        self.fig.get_layout_engine().set(h_pad=0.04, hspace=0.0, w_pad=0.04)
        grid = self.fig.add_gridspec(3, 1, height_ratios=[1.00, 1.00, 1.00])
        self.ax_raw = self.fig.add_subplot(grid[0])
        self.ax = self.fig.add_subplot(grid[1])  # zoomable; the raw plot stays an overview
        self.ax_act = self.fig.add_subplot(grid[2])
        self.canvas = FigureCanvasTkAgg(self.fig, master=centre)
        review = self

        class _Toolbar(NavigationToolbar2Tk):
            def home(self, *args):  # Home button: full view of the detrended plot
                review._full_view()

        self.toolbar = _Toolbar(self.canvas, centre, pack_toolbar=False)
        self.toolbar.pack(side="bottom", fill="x")
        ttk.Label(centre, textvariable=self.reg_var, justify="left", padding=(2, 4),
                  font="TkFixedFont").pack(side="bottom", fill="x")
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self.canvas.mpl_connect("button_press_event", self._on_act_press)
        self.canvas.mpl_connect("motion_notify_event", self._on_act_motion)
        self.canvas.mpl_connect("button_release_event", self._on_act_release)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("button_press_event", self._on_zoom_press)
        self.canvas.mpl_connect("motion_notify_event", self._on_zoom_motion)
        self.canvas.mpl_connect("button_release_event", self._on_zoom_release)

        bottom = ttk.Frame(win, padding=10)
        bottom.grid(row=2, column=0, columnspan=2, sticky="ew")
        ttk.Label(bottom, textvariable=self.status_var).pack(side="left")
        ttk.Button(bottom, text="OK \u2013 export files", command=self._ok).pack(side="right")
        ttk.Button(bottom, text="Cancel", command=self._cancel).pack(side="right", padx=6)
        ttk.Button(bottom, text="Reset all", command=self._reset_all).pack(side="right", padx=6)
        ttk.Button(bottom, text="Load edits...", command=self._load_edits).pack(side="right", padx=6)
        ttk.Button(bottom, text="Save edits...", command=self._save_edits).pack(side="right")

        self._link_to_parent()
        win.bind("<Left>", lambda _e: self._key_step(-1))
        win.bind("<Right>", lambda _e: self._key_step(1))
        win.bind("<comma>", lambda _e: self._key_zoom(lambda: self._jump_marker(-1)))
        win.bind("<period>", lambda _e: self._key_zoom(lambda: self._jump_marker(1)))
        win.bind("<Escape>", lambda _e: self._key_zoom(self._full_view))

        self._draw_full()
        self.status_var.set(f"Well {self.well}.   Click 'Help' for how to edit peaks/troughs "
                            "and choose the regression points.")
        try:
            win.wait_visibility()
            win.grab_set()
        except tk.TclError:
            pass
        win.focus_set()

    # ------------------------------------------------------------------ #
    # Window management (minimize / maximize together with the main window)
    # ------------------------------------------------------------------ #

    def _link_to_parent(self) -> None:
        """Minimizing the review window minimizes the main window too; restoring either
        one brings both back with the review window in front.  While the review window
        is open the main window is blocked, so a click on it raises the review window."""
        self.win.bind("<Unmap>", self._on_review_unmap, add="+")
        self.win.bind("<Map>", self._on_review_map, add="+")
        self.parent.bind("<Map>", self._on_parent_map, add="+")
        self.parent.bind("<FocusIn>", self._on_parent_focus, add="+")

    def _unlink_from_parent(self) -> None:
        for sequence in ("<Map>", "<FocusIn>"):
            try:
                self.parent.unbind(sequence)  # the main window binds nothing else to these
            except tk.TclError:
                pass

    def _window_state(self, window) -> str:
        try:
            return str(window.state())
        except tk.TclError:
            return "withdrawn"

    def _on_review_unmap(self, event) -> None:
        if event.widget is not self.win or self._closed or self._syncing_state:
            return
        if self._window_state(self.win) == "iconic":
            self._syncing_state = True
            try:
                self.parent.iconify()
            except tk.TclError:
                pass
            finally:
                self._syncing_state = False

    def _on_review_map(self, event) -> None:
        if event.widget is not self.win or self._closed or self._syncing_state:
            return
        if self._window_state(self.parent) == "iconic":
            self._syncing_state = True
            try:
                self.parent.deiconify()
            except tk.TclError:
                pass
            finally:
                self._syncing_state = False
        self._raise()

    def _on_parent_map(self, event) -> None:
        if event.widget is not self.parent or self._closed or self._syncing_state:
            return
        if self._window_state(self.win) == "iconic":
            self._syncing_state = True
            try:
                self.win.deiconify()
            except tk.TclError:
                pass
            finally:
                self._syncing_state = False
        self._raise()

    def _on_parent_focus(self, _event=None) -> None:
        if not self._closed and self._window_state(self.win) not in ("iconic", "withdrawn"):
            self.win.after_idle(self._raise)

    def _raise(self) -> None:
        if self._closed:
            return
        try:
            self.win.lift()
            self.win.focus_force()
        except tk.TclError:
            pass

    def _show_help(self) -> None:
        messagebox.showinfo("Review window \u2013 help", self.HELP_TEXT, parent=self.win)

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #

    @property
    def day_length(self) -> float:
        return self.periods[self.well]

    @day_length.setter
    def day_length(self, value: float) -> None:
        self.periods[self.well] = value

    def _all_selected(self) -> dict:
        return {w: (set(p), set(t)) for w, (p, t) in self.selection.items()}

    @staticmethod
    def _clamp_period(value: float) -> float:
        return min(max(round(value, 2), ACTOGRAM_MIN_PERIOD_HOURS), ACTOGRAM_MAX_PERIOD_HOURS)

    def _used(self, well: str | None = None) -> tuple[list[int], list[int]]:
        well = well or self.well
        peaks, troughs = self.selection[well]
        sel_peaks, sel_troughs = self.reg_sel[well]
        return ([i for i in peaks if i in sel_peaks], [i for i in troughs if i in sel_troughs])

    def _fit(self, well: str) -> dict:
        if well not in self._fit_cache:
            self._fit_cache[well] = fit_one_well(self.data, self.analysis.trend.data, well,
                                                 self.analysis.settings)
        return self._fit_cache[well]

    # ------------------------------------------------------------------ #
    # Well list
    # ------------------------------------------------------------------ #

    def _is_edited(self, well: str) -> bool:
        peaks, troughs = self.selection[well]
        auto_peaks, auto_troughs = self.auto[well]
        return sorted(peaks) != sorted(auto_peaks) or sorted(troughs) != sorted(auto_troughs)

    def _list_text(self, well: str) -> str:
        peaks, troughs = self.selection[well]
        mark = "  *" if self._is_edited(well) else ""
        return f"{well:<6} P {len(peaks):>2}  T {len(troughs):>2}{mark}"

    def _refresh_list_item(self, well: str) -> None:
        i = self.wells.index(well)
        self.listbox.delete(i)
        self.listbox.insert(i, self._list_text(well))
        if well == self.well:
            self.listbox.selection_set(i)

    def _on_list_select(self, _event=None) -> None:
        selected = self.listbox.curselection()
        if selected and self.wells[selected[0]] != self.well:
            self._apply_period()
            self.well = self.wells[selected[0]]
            self._draw_full()

    def _key_step(self, delta: int) -> None:
        if self.win.focus_get() is self.period_spin:
            return
        self._step(delta)

    def _step(self, delta: int) -> None:
        self._apply_period()
        i = min(max(self.wells.index(self.well) + delta, 0), len(self.wells) - 1)
        self.listbox.selection_clear(0, "end")
        self.listbox.selection_set(i)
        self.listbox.see(i)
        if self.wells[i] != self.well:
            self.well = self.wells[i]
            self._draw_full()

    # ------------------------------------------------------------------ #
    # Drawing
    # ------------------------------------------------------------------ #

    def _draw_full(self, keep_view: bool = False) -> None:
        """Redraw the well.  keep_view=True keeps the current zoom (e.g. 'Show sine fit')."""
        well = self.well
        a = self.analysis
        ax, ax_raw = self.ax, self.ax_raw
        kept = (ax.get_xlim(), ax.get_ylim()) if keep_view and not self._is_full_view() else None
        ax_raw.clear()
        ax.clear()
        self._marker_artists = []

        ax_raw.scatter(a.raw[TIME], a.raw[well], s=2.0, c="tab:purple", linewidths=0,
                       label="Raw data")
        ax_raw.plot(a.trend.data[TIME], a.trend.data[well], "-r", linewidth=1.0,
                    label=f"Trend: {a.trend.label}")
        ax_raw.set_ylabel(a.settings.value_label.split("(")[0].strip() or "Raw", fontsize=8)
        outside_legend(ax_raw, fontsize=7)
        ax_raw.grid(True, linewidth=0.5, color="lightgray", linestyle="--")
        ax_raw.tick_params(labelbottom=False)
        ax_raw.set_navigate(False)  # zoom/pan applies to the detrended plot only
        self._view_span = None      # ax_raw.clear() above discarded the band
        ax_raw.set_title(f"{a.settings.experiment_name}   Well {well}")

        ax.scatter(self.data[TIME], self.data[well], s=2.5, c="violet", linewidths=0,
                   label="Detrended")
        ax.plot(self.smooth[TIME], self.smooth[well], "-b", linewidth=1.0, label=a.smooth_label)
        fit_text = ""
        if self.fit_var.get():
            detail = self._fit(well)
            if detail["fit"] is not None:
                ax.plot(detail["times"], fit_curve(detail["fit"], detail["times"]), "-",
                        color="darkorange", linewidth=1.2, alpha=0.9, label="Sine fit")
            fit_text = "\n" + fit_summary_text(detail["row"])
        self._fit_text = fit_text
        ax.grid(True, linewidth=0.5, color="lightgray", linestyle="--")
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Detrended")
        ax.autoscale_view()
        low, high = ax.get_ylim()
        pad = (high - low) * 0.10
        ax.set_ylim(low - pad, high + pad)
        hours = self.data[TIME].dropna()
        x_min = math.floor(float(hours.min()) / 12) * 12
        x_max = math.ceil(float(hours.max()) / 12) * 12
        ax.set_xlim(x_min, x_max)
        ax_raw.set_xlim(x_min, x_max)
        ax_raw.xaxis.set_major_locator(MultipleLocator(self._tick_step(x_max - x_min)))
        ax.set_autoscale_on(False)
        self._full_xlim = (x_min, x_max)
        self._full_ylim = ax.get_ylim()
        self._pan = None
        if not keep_view:
            self._focus_row = None
        # ax.clear() dropped earlier callbacks; follow every x change (ours or the toolbar's)
        ax.callbacks.connect("xlim_changed", self._on_view_changed)
        if kept is not None:
            ax.set_xlim(*kept[0])
            ax.set_ylim(*kept[1])
        self._on_view_changed(ax)
        self.period_var.set(f"{self.day_length:g}")
        self.toolbar.update()  # forget the previous well's zoom history
        self._draw_markers()
        self._draw_actogram_preview()
        self.status_var.set(f"Well {well}")

    def _draw_markers(self) -> None:
        for artist in self._marker_artists:
            try:
                artist.remove()
            except (NotImplementedError, ValueError):
                pass
        self._marker_artists = []
        ax = self.ax
        well = self.well
        peaks, troughs = self.selection[well]
        for indices, color, name, sign, va in ((peaks, "red", "Peaks", 1, "bottom"),
                                               (troughs, "blue", "Troughs", -1, "top")):
            if not indices:
                continue
            hours = self.smooth.loc[indices, TIME]
            values = self.smooth.loc[indices, well]
            self._marker_artists.append(
                ax.scatter(hours, values, marker="o", s=45, color=color, zorder=5, label=name))
            for hour, value in zip(hours, values):
                # offset in points, so the labels stay next to the markers at any zoom level
                self._marker_artists.append(ax.annotate(
                    f"{hour:.1f}", (hour, value), xytext=(0, sign * 6),
                    textcoords="offset points", fontsize=7, color=color, ha="center", va=va,
                    annotation_clip=True, clip_on=True))
        focus = self._focus_row
        if focus is not None and (focus in peaks or focus in troughs):
            self._marker_artists.append(ax.scatter(
                [self.smooth.loc[focus, TIME]], [self.smooth.loc[focus, well]], marker="o",
                s=220, facecolors="none", edgecolors="darkorange", linewidths=2.0, zorder=6))
        outside_legend(ax, fontsize=7)

    def _draw_actogram_preview(self) -> None:
        a = self.analysis
        peaks, troughs = self.selection[self.well]
        used = self._used()
        regression = well_regressions(self.smooth, used[0], used[1], self.day_length)
        self._drag_patch = None
        self.ax_act.clear()
        draw_actogram(self.ax_act, self.smooth, peaks, troughs, self.day_length,
                      a.recording.first_hour, a.recording.last_hour, a.settings.label_actogram,
                      selected=used, regression=regression, legend_fontsize=7)
        self.ax_act.set_navigate(False)
        self.reg_var.set(regression_summary_text(regression, len(used[0]), len(used[1]))
                         + getattr(self, "_fit_text", ""))
        self.canvas.draw_idle()

    def _after_edit(self, message: str) -> None:
        self._refresh_list_item(self.well)
        self._draw_markers()
        self._draw_actogram_preview()
        self.status_var.set(message)

    # ------------------------------------------------------------------ #
    # Actogram period
    # ------------------------------------------------------------------ #

    def _set_period(self, value: float) -> None:
        self.period_var.set(f"{self._clamp_period(value):g}")
        self._apply_period()

    def _period_from_fit(self) -> None:
        detail = self._fit(self.well)
        if detail["fit"] is None:
            self.status_var.set("No sine fit is available for this well.")
            return
        self._set_period(round(detail["fit"]["Period (h)"], 1))

    def _apply_period(self, _event=None) -> None:
        if self._closed:
            return
        try:
            value = round(float(self.period_var.get()), 2)
            valid = ACTOGRAM_MIN_PERIOD_HOURS <= value <= ACTOGRAM_MAX_PERIOD_HOURS
        except (ValueError, tk.TclError):
            valid = False
        if not valid:
            self.period_var.set(f"{self.day_length:g}")
            self.status_var.set(f"The actogram period must be between {ACTOGRAM_MIN_PERIOD_HOURS:g} "
                                f"and {ACTOGRAM_MAX_PERIOD_HOURS:g} h.")
            return
        self.period_var.set(f"{value:g}")
        if value != self.day_length:
            self.day_length = value
            self._draw_actogram_preview()
            self.status_var.set(f"Actogram period of well {self.well} set to {value:g} h.")

    def _period_to_all(self) -> None:
        self._apply_period()
        value = self.day_length
        self.periods = {w: value for w in self.wells}
        self._draw_actogram_preview()
        self.status_var.set(f"Actogram period {value:g} h applied to all wells.")

    # ------------------------------------------------------------------ #
    # Editing peaks/troughs in the middle plot
    # ------------------------------------------------------------------ #

    def _toolbar_active(self) -> bool:
        mode = getattr(self.toolbar, "mode", "")
        return bool(getattr(mode, "value", mode))

    def _on_click(self, event) -> None:
        if event.inaxes is not self.ax or event.button != 1 or event.xdata is None:
            return
        if self._toolbar_active():
            return
        mode = self.mode_var.get()
        if mode == "remove":
            self._remove_nearest(event.x, event.y)
        else:
            self._add_at(event.xdata, "peak" if mode == "add_peak" else "trough")

    def _remove_nearest(self, pixel_x: float, pixel_y: float) -> None:
        well = self.well
        peaks, troughs = self.selection[well]
        candidates = [("peak", i) for i in peaks] + [("trough", i) for i in troughs]
        if not candidates:
            self.status_var.set("There are no markers in this well.")
            return
        points = np.array([(self.smooth.loc[i, TIME], self.smooth.loc[i, well])
                           for _, i in candidates], dtype=float)
        pixels = self.ax.transData.transform(points)
        distances = np.hypot(pixels[:, 0] - pixel_x, pixels[:, 1] - pixel_y)
        nearest = int(np.argmin(distances))
        if distances[nearest] > self.HIT_RADIUS_PX:
            self.status_var.set("No marker near the click.")
            return
        kind, index = candidates[nearest]
        (peaks if kind == "peak" else troughs).remove(index)
        if index == self._focus_row:
            self._focus_row = None
        self.reg_sel[well][0 if kind == "peak" else 1].discard(index)
        self._after_edit(f"Removed {kind} at {float(self.smooth.loc[index, TIME]):.2f} h.")

    def _add_at(self, x: float, kind: str) -> None:
        """Add a peak/trough at the local max/min of the smoothed curve near x."""
        well = self.well
        column = self.smooth[well]
        valid = self.smooth.loc[column.notna(), TIME]
        if valid.empty:
            self.status_var.set("This well has no data.")
            return
        near = valid[(valid - x).abs() <= 1.0]  # +/- 1 h search window
        if near.empty:
            index = int((valid - x).abs().idxmin())
        else:
            values = column.loc[near.index]
            index = int(values.idxmax() if kind == "peak" else values.idxmin())
        peaks, troughs = self.selection[well]
        target, other = (peaks, troughs) if kind == "peak" else (troughs, peaks)
        hour = float(self.smooth.loc[index, TIME])
        if index in target:
            self.status_var.set(f"There is already a {kind} at {hour:.2f} h.")
            return
        if index in other:
            self.status_var.set(f"{hour:.2f} h is already marked as the opposite type; remove it first.")
            return
        target.append(index)
        target.sort()
        self.reg_sel[well][0 if kind == "peak" else 1].add(index)
        self._after_edit(f"Added {kind} at {hour:.2f} h.")

    def _reset_well(self) -> None:
        peaks, troughs = self.auto[self.well]
        self.selection[self.well] = (list(peaks), list(troughs))
        self.reg_sel[self.well] = (set(peaks), set(troughs))
        self._after_edit(f"Well {self.well} reset to the automatic detection.")

    def _reset_all(self) -> None:
        if not messagebox.askyesno(APP_TITLE, "Discard all edits in every well?", parent=self.win):
            return
        self.selection = {w: (list(p), list(t)) for w, (p, t) in self.auto.items()}
        self.reg_sel = self._all_selected()
        for well in self.wells:
            self._refresh_list_item(well)
        self._draw_markers()
        self._draw_actogram_preview()
        self.status_var.set("All wells reset to the automatic detection.")

    # ------------------------------------------------------------------ #
    # Zooming the detrended plot with the mouse
    # ------------------------------------------------------------------ #

    def _key_zoom(self, action) -> None:
        # Typing in the period box must not move the view.
        if self.win.focus_get() is self.period_spin:
            return
        action()

    def _is_full_view(self) -> bool:
        x0, x1 = self.ax.get_xlim()
        full0, full1 = self._full_xlim
        return (x1 - x0) >= (full1 - full0) * 0.999

    @staticmethod
    def _tick_step(span: float) -> float:
        """Major tick spacing: 24 h in a long view, finer when zoomed in."""
        for step in (HOURS_PER_DAY, 12.0, 6.0, 4.0, 2.0, 1.0):
            if span / step >= 4:
                return step
        return 0.5

    def _on_view_changed(self, ax) -> None:
        """x range of the detrended plot changed: adapt the ticks, move the band in the raw plot."""
        x0, x1 = ax.get_xlim()
        ax.xaxis.set_major_locator(MultipleLocator(self._tick_step(x1 - x0)))
        if self._view_span is not None:
            try:
                self._view_span.remove()
            except (NotImplementedError, ValueError):
                pass
            self._view_span = None
        if not self._is_full_view():
            self._view_span = self.ax_raw.axvspan(
                x0, x1, facecolor="gold", edgecolor="darkorange", alpha=0.3, zorder=0.5)

    def _fit_y(self, x0: float, x1: float) -> None:
        """Fit the y range to the data inside [x0, x1] (room left for the labels)."""
        well = self.well
        values = []
        for frame in (self.data, self.smooth):
            inside = (frame[TIME] >= x0) & (frame[TIME] <= x1)
            column = frame.loc[inside, well].to_numpy(dtype=float)
            values.append(column[np.isfinite(column)])
        values = np.concatenate(values)
        if values.size == 0:
            return
        low, high = float(values.min()), float(values.max())
        span = high - low
        if span <= 0:
            span = abs(high) or 1.0
        self.ax.set_ylim(low - span * 0.15, high + span * 0.15)

    def _set_view(self, x0: float, x1: float, fit_y: bool = True) -> None:
        """Show [x0, x1] in the detrended plot, kept inside the recording."""
        full0, full1 = self._full_xlim
        full_width = full1 - full0
        min_width = min(full_width,
                        max(self.ZOOM_MIN_WIDTH_HOURS, 10 * self.analysis.time_interval))
        width = min(max(x1 - x0, min_width), full_width)
        if width >= full_width * 0.999:
            x0, x1 = full0, full1
        else:
            x0 = min(max(x0, full0), full1 - width)
            x1 = x0 + width
        self.ax.set_xlim(x0, x1)
        if x0 == full0 and x1 == full1:
            self.ax.set_ylim(*self._full_ylim)
        elif fit_y:
            self._fit_y(x0, x1)
        self.canvas.draw_idle()

    def _zoom_width(self) -> float:
        """Width for zooming onto a point: the current one if already zoomed further in."""
        x0, x1 = self.ax.get_xlim()
        return min(x1 - x0, self.ZOOM_WINDOW_HOURS)

    def _centre_on(self, hour: float, focus_row: int | None = None) -> None:
        width = self._zoom_width()
        if focus_row != self._focus_row:
            self._focus_row = focus_row
            self._draw_markers()
        self._set_view(hour - width / 2, hour + width / 2)

    def _full_view(self) -> None:
        self._focus_row = None
        self._draw_markers()
        self._set_view(*self._full_xlim)
        self.status_var.set(f"Well {self.well}: full view.")

    def _sorted_markers(self) -> list:
        """All peaks and troughs of the current well in time order: (hour, kind, row)."""
        peaks, troughs = self.selection[self.well]
        markers = [(float(self.smooth.loc[i, TIME]), "peak", i) for i in peaks]
        markers += [(float(self.smooth.loc[i, TIME]), "trough", i) for i in troughs]
        return sorted(markers)

    def _jump_marker(self, delta: int) -> None:
        """Zoom onto the next (delta=1) or previous (delta=-1) peak/trough."""
        markers = self._sorted_markers()
        if not markers:
            self.status_var.set("There are no peaks or troughs in this well.")
            return
        x0, x1 = self.ax.get_xlim()
        rows = [row for _h, _k, row in markers]
        if self._is_full_view():
            i = 0 if delta > 0 else len(markers) - 1
        elif self._focus_row in rows and x0 <= self.smooth.loc[self._focus_row, TIME] <= x1:
            i = rows.index(self._focus_row) + delta
        else:  # continue from the middle of the current view
            centre = (x0 + x1) / 2
            if delta > 0:
                i = next((k for k, m in enumerate(markers) if m[0] > centre + 1e-9), len(markers))
            else:
                i = next((k for k in range(len(markers) - 1, -1, -1)
                          if markers[k][0] < centre - 1e-9), -1)
        if not 0 <= i < len(markers):
            self.status_var.set("This is the last peak/trough." if delta > 0
                                else "This is the first peak/trough.")
            return
        hour, kind, row = markers[i]
        self._centre_on(hour, row)
        self.status_var.set(f"{kind.capitalize()} at {hour:.2f} h  ({i + 1} of {len(markers)})")

    def _on_scroll(self, event) -> None:
        if event.inaxes is not self.ax or event.xdata is None:
            return
        state = getattr(getattr(event, "guiEvent", None), "state", 0)
        shift = (isinstance(state, int) and state & 0x0001) or "shift" in (event.key or "")
        x0, x1 = self.ax.get_xlim()
        width = x1 - x0
        if shift:  # scroll in time
            shift_by = -event.step * self.WHEEL_PAN_FRACTION * width
            self._set_view(x0 + shift_by, x1 + shift_by)
            return
        scale = self.WHEEL_ZOOM_FACTOR ** event.step
        anchor = event.xdata  # the time under the pointer stays where it is
        self._set_view(anchor - (anchor - x0) * scale, anchor + (x1 - anchor) * scale)

    def _nearest_marker(self, pixel_x: float, pixel_y: float, radius: float):
        markers = self._sorted_markers()
        if not markers:
            return None
        points = np.array([(hour, self.smooth.loc[row, self.well]) for hour, _k, row in markers],
                          dtype=float)
        pixels = self.ax.transData.transform(points)
        distances = np.hypot(pixels[:, 0] - pixel_x, pixels[:, 1] - pixel_y)
        nearest = int(np.argmin(distances))
        return markers[nearest] if distances[nearest] <= radius else None

    def _on_zoom_press(self, event) -> None:
        if self._toolbar_active() or event.xdata is None:
            return
        if event.inaxes is self.ax_raw and event.button == 1:
            # the raw-data plot works as an overview: jump there
            self._centre_on(event.xdata)
            self.status_var.set(f"View moved to {event.xdata:.1f} h.")
            return
        if event.inaxes is not self.ax:
            return
        if event.button == 2 or (event.button == 3 and getattr(event, "dblclick", False)):
            self._pan = None
            self._full_view()
        elif event.button == 3:
            self._pan = {"px": event.x, "py": event.y, "xlim": self.ax.get_xlim(), "moved": False}

    def _on_zoom_motion(self, event) -> None:
        pan = self._pan
        if pan is None or event.x is None:
            return
        if not pan["moved"] and abs(event.x - pan["px"]) < self.DRAG_THRESHOLD_PX:
            return
        pan["moved"] = True
        x0, x1 = pan["xlim"]
        hours_per_pixel = (x1 - x0) / max(self.ax.bbox.width, 1.0)
        shift_by = -(event.x - pan["px"]) * hours_per_pixel
        self._set_view(x0 + shift_by, x1 + shift_by, fit_y=False)

    def _on_zoom_release(self, event) -> None:
        pan, self._pan = self._pan, None
        if pan is None:
            return
        if pan["moved"]:
            x0, x1 = self.ax.get_xlim()
            self._set_view(x0, x1)  # fit the y range to what is now visible
            return
        # a right-click without dragging: zoom onto the nearest peak/trough
        marker = self._nearest_marker(pan["px"], pan["py"], self.HIT_RADIUS_PX * 3)
        if marker is None:
            x = self.ax.transData.inverted().transform((pan["px"], 0))[0]
            self._centre_on(x)
            self.status_var.set(f"Zoomed in at {x:.1f} h.")
            return
        hour, kind, row = marker
        self._centre_on(hour, row)
        markers = self._sorted_markers()
        self.status_var.set(f"{kind.capitalize()} at {hour:.2f} h  "
                            f"({markers.index(marker) + 1} of {len(markers)})")

    # ------------------------------------------------------------------ #
    # Choosing the regression points in the actogram
    # ------------------------------------------------------------------ #

    def _actogram_hits(self) -> list:
        peaks, troughs = self.selection[self.well]
        hits = []
        for kind_number, indices in enumerate((peaks, troughs)):
            for x, y, _hour, row in actogram_points(self.smooth, indices, self.day_length):
                hits.append((kind_number, row, x, y))
        return hits

    def _set_selected(self, kind_number: int, row: int, selected: bool) -> None:
        target = self.reg_sel[self.well][kind_number]
        if selected:
            target.add(row)
        else:
            target.discard(row)

    def _select_all(self, selected: bool) -> None:
        peaks, troughs = self.selection[self.well]
        self.reg_sel[self.well] = (set(peaks), set(troughs)) if selected else (set(), set())
        self._draw_actogram_preview()
        self.status_var.set("All peaks and troughs of this well are used in the regression."
                            if selected else "No points selected for the regression in this well.")

    def _on_act_press(self, event) -> None:
        if event.inaxes is not self.ax_act or event.button not in (1, 3):
            return
        if event.xdata is None or self._toolbar_active():
            return
        self._drag = {"x": event.xdata, "y": event.ydata, "px": event.x, "py": event.y,
                      "button": event.button}

    def _on_act_motion(self, event) -> None:
        drag = self._drag
        if drag is None or event.x is None:
            return
        if math.hypot(event.x - drag["px"], event.y - drag["py"]) < self.DRAG_THRESHOLD_PX:
            return
        x1, y1 = self.ax_act.transData.inverted().transform((event.x, event.y))
        x0, y0 = drag["x"], drag["y"]
        if self._drag_patch is None:
            from matplotlib.patches import Rectangle

            color = "tab:green" if drag["button"] == 1 else "tab:gray"
            self._drag_patch = self.ax_act.add_patch(Rectangle(
                (x0, y0), 0, 0, fill=True, alpha=0.2, facecolor=color, edgecolor=color,
                linestyle="--", zorder=6))
        self._drag_patch.set_bounds(min(x0, x1), min(y0, y1), abs(x1 - x0), abs(y1 - y0))
        self.canvas.draw_idle()

    def _on_act_release(self, event) -> None:
        drag, self._drag = self._drag, None
        if self._drag_patch is not None:
            try:
                self._drag_patch.remove()
            except (NotImplementedError, ValueError):
                pass
            self._drag_patch = None
        if drag is None or event.x is None:
            return
        moved = math.hypot(event.x - drag["px"], event.y - drag["py"])
        if moved < self.DRAG_THRESHOLD_PX:
            if drag["button"] == 1:
                self._toggle_nearest(drag["px"], drag["py"])
            else:
                self.canvas.draw_idle()
            return
        x1, y1 = self.ax_act.transData.inverted().transform((event.x, event.y))
        x_lo, x_hi = sorted((drag["x"], x1))
        y_lo, y_hi = sorted((drag["y"], y1))
        select = drag["button"] == 1
        count = 0
        for kind_number, row, x, y in self._actogram_hits():
            if x_lo <= x <= x_hi and y_lo <= y <= y_hi:
                self._set_selected(kind_number, row, select)
                count += 1
        self._draw_actogram_preview()
        self.status_var.set(f"{'Selected' if select else 'Deselected'} the points in the box "
                            f"({count} drawn points touched).")

    def _toggle_nearest(self, pixel_x: float, pixel_y: float) -> None:
        hits = self._actogram_hits()
        if not hits:
            self.status_var.set("There are no peaks or troughs in this well.")
            return
        points = np.array([(x, y) for _k, _r, x, y in hits], dtype=float)
        pixels = self.ax_act.transData.transform(points)
        distances = np.hypot(pixels[:, 0] - pixel_x, pixels[:, 1] - pixel_y)
        nearest = int(np.argmin(distances))
        if distances[nearest] > self.HIT_RADIUS_PX:
            self.status_var.set("No point near the click.")
            return
        kind_number, row, _x, _y = hits[nearest]
        now_selected = row not in self.reg_sel[self.well][kind_number]
        self._set_selected(kind_number, row, now_selected)
        self._draw_actogram_preview()
        hour = float(self.smooth.loc[row, TIME])
        kind = "peak" if kind_number == 0 else "trough"
        self.status_var.set(f"{'Using' if now_selected else 'Not using'} the {kind} at "
                            f"{hour:.2f} h for the regression.")

    # ------------------------------------------------------------------ #
    # Save / load / finish
    # ------------------------------------------------------------------ #

    def _save_edits(self) -> None:
        name = safe_filename(self.analysis.settings.experiment_name) or "experiment"
        path = filedialog.asksaveasfilename(
            parent=self.win, title="Save peak/trough edits", defaultextension=".json",
            initialfile=f"{name}_peaks.json", filetypes=[("JSON files", "*.json")])
        if not path:
            return
        self._apply_period()
        try:
            save_peaks_json(path, self.selection, self.smooth[TIME], self.analysis.input_path.name,
                            {w: self._used(w) for w in self.wells}, self.periods)
        except OSError as error:
            messagebox.showerror(APP_TITLE, f"Could not save the file: {error}", parent=self.win)
            return
        self.status_var.set(f"Saved: {Path(path).name}")

    def _load_edits(self) -> None:
        path = filedialog.askopenfilename(
            parent=self.win, title="Load peak/trough edits",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if not path:
            return
        try:
            loaded, skipped, loaded_reg, loaded_period = load_peaks_json(
                path, self.smooth[TIME], self.wells)
        except (OSError, ValueError) as error:
            messagebox.showerror(APP_TITLE, f"Could not load the file: {error}", parent=self.win)
            return
        self.selection.update(loaded)
        for well, (peaks, troughs) in loaded.items():
            self.reg_sel[well] = (set(peaks), set(troughs))
        for well, (peaks, troughs) in (loaded_reg or {}).items():
            if well in loaded:
                self.reg_sel[well] = (set(peaks), set(troughs))
        for well, value in (loaded_period or {}).items():
            self.periods[well] = self._clamp_period(value)
        for well in self.wells:
            self._refresh_list_item(well)
        self._draw_full()
        note = (f" ({skipped} time points outside this analysis were skipped)" if skipped else "")
        self.status_var.set(f"Loaded {len(loaded)} wells from {Path(path).name}{note}")

    def _close(self) -> None:
        self._closed = True
        self._unlink_from_parent()
        try:
            self.win.grab_release()
        except tk.TclError:
            pass
        self.win.destroy()

    def _ok(self) -> None:
        self._apply_period()
        selection = {w: (sorted(p), sorted(t)) for w, (p, t) in self.selection.items()}
        reg_selection = {w: self._used(w) for w in self.wells}
        periods = dict(self.periods)
        self._close()
        self.on_ok(selection, reg_selection, periods)

    def _cancel(self) -> None:
        self._close()
        self.on_cancel()


def _ensure_std_streams() -> None:
    """A windowed EXE (PyInstaller --windowed) has no console: sys.stdout/stderr are
    None, and libraries that write to them would fail.  Send such output nowhere."""
    for name in ("stdout", "stderr"):
        if getattr(sys, name) is None:
            setattr(sys, name, open(os.devnull, "w", encoding="utf-8"))


def self_test(out_dir: str = "selftest") -> int:
    """Headless check that every part works (used after building the EXE).

    Makes a small synthetic workbook (12 wells, known periods), runs the whole
    analysis including the Excel/PDF output, checks the results and writes
    '<out_dir>/selftest.log'.  Returns the process exit code (0 = OK).
    """
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "selftest.log"
    with open(log_path, "w", encoding="utf-8") as log, \
            contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        try:
            print(f"{APP_TITLE} {APP_VERSION} self-test, {utc_now_string()} UTC")
            print(f"Python {sys.version.split()[0]}, numpy {np.__version__}, "
                  f"pandas {pd.__version__}, matplotlib {matplotlib.__version__}, "
                  f"frozen={getattr(sys, 'frozen', False)}")
            if tk is None:
                raise RuntimeError("tkinter is not available")
            from matplotlib.backends import backend_tkagg  # noqa: F401 - bundled?

            rng = np.random.default_rng(0)
            hours = np.arange(0, 96.01, 0.25)
            truth = {}
            columns = {"Time (h)": hours, "Temperature (\u00b0C)": 35 + rng.normal(0, 0.05, hours.size)}
            for i, well in enumerate(f"{r}{c}" for r in "CD" for c in range(1, 7)):
                period = 20.0 + 0.5 * i
                truth[well] = period
                columns[well] = (110 - 0.05 * hours
                                 + 6 * np.exp(-hours / 60) * np.cos(2 * np.pi * hours / period)
                                 + rng.normal(0, 0.5, hours.size))
            source = out / "selftest_input.xlsx"
            pd.DataFrame(columns).to_excel(source, index=False)

            settings = Settings(experiment_name="selftest", fit_start_hour=0.0)
            files = run_analysis(source, out, settings)
            missing = [f.name for f in files if not Path(f).is_file()]
            if missing:
                raise RuntimeError(f"output files missing: {missing}")
            summary = pd.read_excel(files[0], sheet_name="Summary").set_index("Well")
            errors = {w: abs(summary.loc[w, "Fit Period (h)"] - p) for w, p in truth.items()}
            worst = max(errors, key=errors.get)
            print(f"Largest fitted-period error: {errors[worst]:.3f} h (well {worst})")
            if errors[worst] > 0.5:
                raise RuntimeError("fitted periods do not match the synthetic data")
            print("SELFTEST OK")
            return 0
        except Exception:  # noqa: BLE001 - everything goes to the log
            traceback.print_exc()
            print("SELFTEST FAILED")
            return 1


def main() -> None:
    _ensure_std_streams()
    args = sys.argv[1:]
    if args and args[0] == "--version":
        print(f"{APP_TITLE} {APP_VERSION}")
        return
    if args and args[0] == "--selftest":
        sys.exit(self_test(args[1] if len(args) > 1 else "selftest"))
    if tk is None:
        sys.exit("Tkinter is not available in this Python installation.")
    if sys.platform.startswith("win"):
        try:
            import ctypes

            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:  # noqa: BLE001
            pass
    root = tk.Tk()
    if sys.platform.startswith("linux"):
        style = ttk.Style(root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
