"""Generic scatter helpers: per-model, RSCC residual, aggregated."""

from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text

from style import MODEL_DISPLAY, MODEL_PALETTE


def make_scatter_comparison(
    df: pd.DataFrame,
    model_name: str,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    title: str | None = None,
    delta_col: str | None = None,
    draw_diagonal: bool = False,
    highlight_top_n: int = 3,
    highlight_decreased: bool | int = True,
    highlight_proteins: list[str] | None = None,
    label_col: str = "protein",
    label_all: bool = False,
    draw_regression: bool = False,
    show_correlation: bool = False,
    equal_aspect: bool = True,
    show_marginal_kde: bool = False,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    residual: bool = False,
    label_top_n: int = 3,
    label_bottom_n: int = 5,
    label_far_from_ceiling_n: int = 5,
) -> plt.Figure:
    """Generic scatter plot comparing two metrics for a single model.

    See module docstring of ``analysis_helpers`` (pre-refactor) or commit
    history for the full prose; behavior is unchanged.
    """
    color = MODEL_PALETTE.get(model_name, "#888888")
    display_name = MODEL_DISPLAY.get(model_name, model_name)
    if title is None:
        title = f"{display_name}: {y_label} vs {x_label}"

    _w = 8 if label_all else 6
    if show_marginal_kde:
        fig = plt.figure(figsize=(_w, _w * 1.25))
        gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
        ax = fig.add_subplot(gs[0])
        ax_kde = fig.add_subplot(gs[1], sharex=ax)
    else:
        fig, ax = plt.subplots(figsize=(_w, _w))
        ax_kde = None

    x_raw = df[x_col]
    y_raw = df[y_col]
    if residual:
        y_plot = y_raw - x_raw
        if "Improvement" not in y_label:
            y_label_used = "Improvement (Guided − Baseline RSCC)"
        else:
            y_label_used = y_label
        _delta_for_highlight = y_plot
        equal_aspect = False
        draw_diagonal_used = False
    else:
        y_plot = y_raw
        y_label_used = y_label
        _delta_for_highlight = df[delta_col] if (delta_col is not None and delta_col in df.columns) else None
        draw_diagonal_used = draw_diagonal

    x = x_raw

    if draw_diagonal_used:
        lo = min(x.min(), y_plot.min()) - 0.02
        hi = max(x.max(), y_plot.max()) + 0.02
        if xlim is not None:
            lo = min(lo, xlim[0])
            hi = max(hi, xlim[1])
        if ylim is not None:
            lo = min(lo, ylim[0])
            hi = max(hi, ylim[1])
        ax.plot([lo, hi], [lo, hi], ls="--", color="#AAAAAA", lw=1, zorder=1)

    if residual:
        _x_lo = float(np.nanmin(x)) if xlim is None else xlim[0]
        _x_hi = float(np.nanmax(x)) if xlim is None else xlim[1]
        _x_lo = min(_x_lo, 0.0) - 0.02
        _x_hi = max(_x_hi, 1.0) + 0.02
        ax.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=1)
        _xs = np.linspace(_x_lo, _x_hi, 100)
        ax.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=1)

    if residual:
        ax.scatter(x, y_plot, c=color, s=14, alpha=0.3, edgecolors="none", zorder=2)
    else:
        ax.scatter(x, y_plot, c=color, s=40, alpha=0.7, edgecolors="white", linewidths=0.5, zorder=2)

    texts: list = []

    if residual:
        _has_label_col = label_col in df.columns
        _seen_idx: set = set()

        if label_top_n > 0 and _delta_for_highlight is not None:
            imp = df.assign(_d=_delta_for_highlight).nlargest(label_top_n, "_d")
            imp = imp[imp["_d"] > 0]
            if not imp.empty:
                ax.scatter(
                    imp[x_col], imp[y_col] - imp[x_col],
                    facecolors=color, s=24, edgecolors="#2E7D32", linewidths=1.2, zorder=3,
                )
                if _has_label_col:
                    for _, r in imp.iterrows():
                        _seen_idx.add(r.name)
                        texts.append(ax.text(
                            r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if _delta_for_highlight is not None:
            dec = df.assign(_d=_delta_for_highlight)
            dec = dec[dec["_d"] < 0]
            dec = dec.loc[~dec.index.isin(_seen_idx)]
            if not dec.empty:
                ax.scatter(
                    dec[x_col], dec[y_col] - dec[x_col],
                    facecolors=color, s=24, edgecolors="#D32F2F", linewidths=1.2, zorder=3,
                )
                if _has_label_col:
                    for _, r in dec.iterrows():
                        _seen_idx.add(r.name)
                        texts.append(ax.text(
                            r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if label_far_from_ceiling_n > 0:
            below_dist = 1.0 - df[y_col]
            far = df.assign(_b=below_dist).nlargest(label_far_from_ceiling_n, "_b")
            far = far[far["_b"] > 0]
            far = far.loc[~far.index.isin(_seen_idx)]
            if not far.empty:
                ax.scatter(
                    far[x_col], far[y_col] - far[x_col],
                    facecolors=color, s=24, edgecolors="#F57C00", linewidths=1.2, zorder=3,
                )
                if _has_label_col:
                    for _, r in far.iterrows():
                        _seen_idx.add(r.name)
                        texts.append(ax.text(
                            r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if highlight_proteins and _has_label_col:
            wanted = {p.lower() for p in highlight_proteins}
            picked = df[df[label_col].str.lower().isin(wanted)]
            if not picked.empty:
                ax.scatter(
                    picked[x_col], picked[y_col] - picked[x_col],
                    facecolors=color, s=24, edgecolors="#1565C0", linewidths=1.2, zorder=4,
                )
                for _, r in picked.iterrows():
                    if r.name in _seen_idx:
                        continue
                    _seen_idx.add(r.name)
                    texts.append(ax.text(
                        r[x_col], r[y_col] - r[x_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

    else:
        if delta_col is not None and delta_col in df.columns:
            if highlight_decreased:
                dec = df[df[delta_col] < 0]
                if isinstance(highlight_decreased, bool):
                    pass
                else:
                    dec = dec.nsmallest(int(highlight_decreased), delta_col)
                if not dec.empty:
                    ax.scatter(
                        dec[x_col], dec[y_col],
                        facecolors=color, s=60, edgecolors="#D32F2F", linewidths=1.5, zorder=3,
                    )
                    for _, r in dec.iterrows():
                        texts.append(ax.text(
                            r[x_col], r[y_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

            if highlight_top_n > 0:
                imp = df.nlargest(highlight_top_n, delta_col)
                imp = imp[imp[delta_col] > 0]
                if not imp.empty:
                    ax.scatter(
                        imp[x_col], imp[y_col],
                        facecolors=color, s=60, edgecolors="#2E7D32", linewidths=1.5, zorder=3,
                    )
                    for _, r in imp.iterrows():
                        texts.append(ax.text(
                            r[x_col], r[y_col], str(r[label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if highlight_proteins and label_col in df.columns:
            wanted = {p.lower() for p in highlight_proteins}
            picked = df[df[label_col].str.lower().isin(wanted)]
            if not picked.empty:
                ax.scatter(
                    picked[x_col], picked[y_col],
                    facecolors=color, s=60, edgecolors="#1565C0", linewidths=1.5, zorder=4,
                )
                for _, r in picked.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        if label_all and delta_col is None:
            for _, r in df.iterrows():
                texts.append(ax.text(
                    r[x_col], r[y_col], str(r[label_col]).upper(),
                    color="#333333", zorder=10,
                ))

    _finite = np.isfinite(x) & np.isfinite(y_plot)
    _n_dropped = int((~_finite).sum())
    if _n_dropped > 0:
        warnings.warn(
            f"make_scatter_comparison({display_name}): dropping {_n_dropped} of "
            f"{len(df)} rows with NaN in '{x_col}' or '{y_col}' before regression.",
            stacklevel=2,
        )
    _xf, _yf = x[_finite], y_plot[_finite]

    if draw_regression and len(_xf) >= 2 and _xf.min() != _xf.max():
        m, b = np.polyfit(_xf, _yf, 1)
        x_range = np.linspace(_xf.min(), _xf.max(), 100)
        ax.plot(x_range, m * x_range + b, ls="-", color=color, lw=1.5, alpha=0.6, zorder=1)

    if show_correlation and len(_xf) >= 5:
        r_val = np.corrcoef(_xf, _yf)[0, 1]
        ax.text(
            0.97, 0.03, f"r = {r_val:.3f}",
            transform=ax.transAxes, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.8),
        )

    ax.set_ylabel(y_label_used)
    ax.set_title(title, pad=14)

    if ax_kde is not None:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis="x", length=0)
        sns.kdeplot(x_raw, ax=ax_kde, color="#888888", fill=True, alpha=0.3,
                    label="Baseline", lw=1.5, bw_adjust=0.5)
        sns.kdeplot(y_raw, ax=ax_kde, color=color, fill=True, alpha=0.3,
                    label="Guided", lw=1.5, bw_adjust=0.5)
        ax_kde.set_xlabel(x_label)
        ax_kde.set_ylabel("Density")
        ax_kde.set_yticks([])
        ax_kde.legend(loc="upper right", frameon=True, framealpha=0.8, edgecolor="none")
    else:
        ax.set_xlabel(x_label)

    _fixed_limits = xlim is not None or ylim is not None
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box" if _fixed_limits else "datalim")
    if equal_aspect and draw_diagonal_used and not _fixed_limits:
        pad = 0.01
        ax.set_xlim(lo + pad, hi)
        ax.set_ylim(lo + pad, hi)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if residual:
        if xlim is None:
            ax.set_xlim(0.0, 1.0)
        if ylim is None:
            _yvals = y_plot[np.isfinite(y_plot)]
            if len(_yvals) > 0:
                _y_lo = float(min(0.0, np.nanmin(_yvals))) - 0.02
            else:
                _y_lo = -0.05
            ax.set_ylim(_y_lo, 1.0)

    if texts:
        try:
            _n = len(texts)
            if _n > 15:
                adjust_text(
                    texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                    expand=(1.5, 1.8),
                    force_text=(1.0, 1.5),
                    force_static=(0.8, 1.2),
                    force_explode=(0.5, 0.8),
                    ensure_inside_axes=False,
                    max_move=(40, 40),
                    time_lim=10,
                )
            else:
                adjust_text(
                    texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                    expand=(1.3, 1.5),
                    force_text=(0.5, 0.8),
                    force_static=(0.5, 0.8),
                    force_explode=(0.3, 0.6),
                    ensure_inside_axes=True,
                    max_move=(20, 20),
                )
        except Exception as exc:
            warnings.warn(
                f"make_scatter_comparison({display_name}): adjust_text failed "
                f"({exc}); dropping labels rather than overlapping.",
                stacklevel=2,
            )
            for t in texts:
                t.remove()

    fig.tight_layout()
    return fig


def make_rscc_scatter(
    df: pd.DataFrame,
    model_name: str,
    optimal_gw: float,
) -> plt.Figure:
    """Per-model residual plot: baseline RSCC vs improvement (guided − baseline).

    ``df`` must contain a ``label`` column formatted ``"<PDB> <range>"``
    as produced by the RSCC scatter cell in ``02_fit_quality.py``. Four
    label sets are drawn: top-3 most improved, every decreased point, the
    5 points furthest below the y = 1 − x ceiling, and the manually
    highlighted ``6B8X 181-183`` and ``5MHX 162-164`` selections.
    """
    return make_scatter_comparison(
        df,
        model_name,
        x_col="rscc_baseline",
        y_col="rscc_optimal",
        x_label="Baseline RSCC",
        y_label="Improvement (Guided − Baseline RSCC)",
        title=f"{MODEL_DISPLAY.get(model_name, model_name)}: RSCC Improvement Over Baseline",
        delta_col="delta_rscc",
        label_col="label",
        show_marginal_kde=True,
        xlim=(0.0, 1.0),
        residual=True,
        label_top_n=3,
        label_far_from_ceiling_n=5,
        highlight_proteins=["6B8X 181-183", "5MHX 162-164"],
    )


def make_aggregated_scatter(
    agg_df: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    x_q1_col: str | None = None,
    x_q3_col: str | None = None,
    y_q1_col: str | None = None,
    y_q3_col: str | None = None,
    x_label: str,
    y_label: str,
    title: str,
    label_col: str = "protein",
    draw_diagonal: bool = False,
    equal_aspect: bool = False,
    point_color: str = "#4D4D4D",
    label_points: bool = False,
    highlight_top_n: int = 3,
    highlight_decreased: bool | int = True,
    highlight_proteins: list[str] | None = None,
    show_correlation: bool = False,
    show_marginal_kde: bool = False,
    show_errors: bool = True,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    residual: bool = False,
    label_top_n: int = 3,
    label_bottom_n: int | None = None,
    label_far_from_ceiling_n: int = 5,
) -> plt.Figure:
    """Scatter plot with asymmetric IQR error bars for cross-model aggregations.

    See pre-refactor ``analysis_helpers`` docstring or commit history for the
    full prose; behavior unchanged.
    """
    if show_marginal_kde:
        fig = plt.figure(figsize=(6, 7.5))
        gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
        ax = fig.add_subplot(gs[0])
        ax_kde = fig.add_subplot(gs[1], sharex=ax)
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax_kde = None

    x_med = agg_df[x_col].to_numpy()
    y_med_raw = agg_df[y_col].to_numpy()

    _err_cols_present = all(
        c is not None for c in (x_q1_col, x_q3_col, y_q1_col, y_q3_col)
    )
    draw_errors = show_errors and _err_cols_present

    if draw_errors:
        x_q1 = agg_df[x_q1_col].to_numpy()
        x_q3 = agg_df[x_q3_col].to_numpy()
        y_q1 = agg_df[y_q1_col].to_numpy()
        y_q3 = agg_df[y_q3_col].to_numpy()

    if residual:
        y_med = y_med_raw - x_med
        equal_aspect = False
        if "Improvement" not in y_label:
            y_label_used = "Improvement (Guided − Baseline RSCC)"
        else:
            y_label_used = y_label
        if draw_errors:
            sd_x = (x_q3 - x_q1) / 2.0
            sd_y = (y_q3 - y_q1) / 2.0
            sd_resid = np.sqrt(sd_y**2 + sd_x**2)
            xerr = np.vstack([x_med - x_q1, x_q3 - x_med])
            yerr = np.vstack([sd_resid, sd_resid])
    else:
        y_med = y_med_raw
        y_label_used = y_label
        if draw_errors:
            xerr = np.vstack([x_med - x_q1, x_q3 - x_med])
            yerr = np.vstack([y_med - y_q1, y_q3 - y_med])

    if draw_diagonal and not residual:
        if draw_errors:
            lo = float(min(np.nanmin(x_q1), np.nanmin(y_q1))) - 0.02
            hi = float(max(np.nanmax(x_q3), np.nanmax(y_q3))) + 0.02
        else:
            lo = float(min(np.nanmin(x_med), np.nanmin(y_med))) - 0.02
            hi = float(max(np.nanmax(x_med), np.nanmax(y_med))) + 0.02
        if xlim is not None:
            lo = min(lo, xlim[0])
            hi = max(hi, xlim[1])
        if ylim is not None:
            lo = min(lo, ylim[0])
            hi = max(hi, ylim[1])
        ax.plot([lo, hi], [lo, hi], ls="--", color="#AAAAAA", lw=1, zorder=1)

    if residual:
        _x_lo = float(np.nanmin(x_med)) if xlim is None else xlim[0]
        _x_hi = float(np.nanmax(x_med)) if xlim is None else xlim[1]
        _x_lo = min(_x_lo, 0.0) - 0.02
        _x_hi = max(_x_hi, 1.0) + 0.02
        ax.axhline(0.0, ls="--", color="#AAAAAA", lw=1, zorder=1)
        _xs = np.linspace(_x_lo, _x_hi, 100)
        ax.plot(_xs, 1.0 - _xs, ls="--", color="#AAAAAA", lw=0.8, zorder=1)

    if draw_errors:
        ax.errorbar(
            x_med, y_med,
            xerr=xerr, yerr=yerr,
            fmt="none",
            ecolor=point_color,
            alpha=0.35,
            capsize=3,
            lw=1.0,
            zorder=2,
        )
    if residual:
        ax.scatter(
            x_med, y_med,
            c=point_color, s=18, alpha=0.3,
            edgecolors="none", zorder=4,
        )
    else:
        ax.scatter(
            x_med, y_med,
            c=point_color, s=55,
            edgecolors="white", linewidths=0.5, zorder=4,
        )

    texts: list = []
    _delta = agg_df[y_col] - agg_df[x_col]

    def _yp(idx):
        if residual:
            return agg_df.loc[idx, y_col] - agg_df.loc[idx, x_col]
        return agg_df.loc[idx, y_col]

    if residual:
        _seen_idx: set = set()

        if label_top_n > 0:
            imp_idx = _delta.nlargest(label_top_n).index
            imp_idx = [i for i in imp_idx if _delta.loc[i] > 0]
            if imp_idx:
                imp = agg_df.loc[imp_idx]
                ax.scatter(
                    imp[x_col], _yp(imp.index),
                    facecolors=point_color, s=32,
                    edgecolors="#2E7D32", linewidths=1.2, zorder=5,
                )
                if label_col in agg_df.columns:
                    for i in imp_idx:
                        _seen_idx.add(i)
                        texts.append(ax.text(
                            agg_df.loc[i, x_col], _yp(i),
                            str(agg_df.loc[i, label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        dec_mask = _delta < 0
        if label_bottom_n is not None and label_bottom_n > 0:
            dec_ranked = _delta[dec_mask].nsmallest(label_bottom_n).index
            dec_idx = [i for i in dec_ranked if i not in _seen_idx]
        else:
            dec_idx = [i for i in agg_df.index[dec_mask] if i not in _seen_idx]
        if dec_idx:
            dec = agg_df.loc[dec_idx]
            ax.scatter(
                dec[x_col], _yp(dec.index),
                facecolors=point_color, s=32,
                edgecolors="#D32F2F", linewidths=1.2, zorder=5,
            )
            if label_col in agg_df.columns:
                for i in dec_idx:
                    _seen_idx.add(i)
                    texts.append(ax.text(
                        agg_df.loc[i, x_col], _yp(i),
                        str(agg_df.loc[i, label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        if label_far_from_ceiling_n > 0:
            below_dist = 1.0 - agg_df[y_col]
            far_idx = below_dist.nlargest(label_far_from_ceiling_n).index
            far_idx = [i for i in far_idx if below_dist.loc[i] > 0 and i not in _seen_idx]
            if far_idx:
                far = agg_df.loc[far_idx]
                ax.scatter(
                    far[x_col], _yp(far.index),
                    facecolors=point_color, s=32,
                    edgecolors="#F57C00", linewidths=1.2, zorder=5,
                )
                if label_col in agg_df.columns:
                    for i in far_idx:
                        _seen_idx.add(i)
                        texts.append(ax.text(
                            agg_df.loc[i, x_col], _yp(i),
                            str(agg_df.loc[i, label_col]).upper(),
                            color="#333333", zorder=10,
                        ))

        if highlight_proteins and label_col in agg_df.columns:
            wanted = {p.lower() for p in highlight_proteins}
            picked_idx = [
                i for i in agg_df.index
                if str(agg_df.loc[i, label_col]).lower() in wanted
            ]
            if picked_idx:
                picked = agg_df.loc[picked_idx]
                ax.scatter(
                    picked[x_col], _yp(picked.index),
                    facecolors=point_color, s=32,
                    edgecolors="#1565C0", linewidths=1.2, zorder=6,
                )
                for i in picked_idx:
                    if i in _seen_idx:
                        continue
                    _seen_idx.add(i)
                    texts.append(ax.text(
                        agg_df.loc[i, x_col], _yp(i),
                        str(agg_df.loc[i, label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

    else:
        if highlight_decreased and label_col in agg_df.columns:
            dec_mask = _delta < 0
            dec = agg_df[dec_mask]
            if not isinstance(highlight_decreased, bool):
                dec_delta = _delta[dec_mask]
                dec = dec.loc[dec_delta.nsmallest(int(highlight_decreased)).index]
            if not dec.empty:
                ax.scatter(
                    dec[x_col], dec[y_col],
                    facecolors=point_color, s=70,
                    edgecolors="#D32F2F", linewidths=1.5, zorder=5,
                )
                for _, r in dec.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        if highlight_top_n > 0 and label_col in agg_df.columns:
            imp = agg_df.loc[_delta.nlargest(highlight_top_n).index]
            imp = imp[_delta.loc[imp.index] > 0]
            if not imp.empty:
                ax.scatter(
                    imp[x_col], imp[y_col],
                    facecolors=point_color, s=70,
                    edgecolors="#2E7D32", linewidths=1.5, zorder=5,
                )
                for _, r in imp.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

        if highlight_proteins and label_col in agg_df.columns:
            wanted = {p.lower() for p in highlight_proteins}
            picked = agg_df[agg_df[label_col].str.lower().isin(wanted)]
            already = {t.get_text().lower() for t in texts}
            picked = picked[~picked[label_col].str.upper().str.lower().isin(already)]
            if not picked.empty:
                ax.scatter(
                    picked[x_col], picked[y_col],
                    facecolors=point_color, s=70,
                    edgecolors="#1565C0", linewidths=1.5, zorder=5,
                )
                for _, r in picked.iterrows():
                    texts.append(ax.text(
                        r[x_col], r[y_col], str(r[label_col]).upper(),
                        color="#333333", zorder=10,
                    ))

    if show_correlation:
        _finite = np.isfinite(x_med) & np.isfinite(y_med)
        if int(_finite.sum()) >= 5:
            r_val = float(np.corrcoef(x_med[_finite], y_med[_finite])[0, 1])
            ax.text(
                0.97, 0.03, f"r = {r_val:.3f}",
                transform=ax.transAxes, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.8),
            )

    ax.set_ylabel(y_label_used)
    ax.set_title(title, pad=14)

    if ax_kde is not None:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis="x", length=0)
        sns.kdeplot(x_med, ax=ax_kde, color="#888888", fill=True, alpha=0.3,
                    label="Baseline", lw=1.5, bw_adjust=0.5)
        sns.kdeplot(y_med_raw, ax=ax_kde, color=point_color, fill=True, alpha=0.3,
                    label="Guided", lw=1.5, bw_adjust=0.5)
        ax_kde.set_xlabel(x_label)
        ax_kde.set_ylabel("Density")
        ax_kde.set_yticks([])
        ax_kde.legend(loc="upper right", frameon=True, framealpha=0.8, edgecolor="none")
    else:
        ax.set_xlabel(x_label)

    _fixed_limits = xlim is not None or ylim is not None
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box" if _fixed_limits else "datalim")
    if equal_aspect and draw_diagonal and not residual and not _fixed_limits:
        pad = 0.01
        ax.set_xlim(lo + pad, hi)
        ax.set_ylim(lo + pad, hi)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if residual:
        if xlim is None:
            ax.set_xlim(0.0, 1.0)
        if ylim is None:
            _yvals = y_med[np.isfinite(y_med)]
            if len(_yvals) > 0:
                _y_lo = float(min(0.0, np.nanmin(_yvals))) - 0.02
            else:
                _y_lo = -0.05
            ax.set_ylim(_y_lo, 1.0)

    if texts:
        try:
            adjust_text(
                texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                expand=(1.3, 1.5),
                force_text=(0.5, 0.8),
                force_static=(0.5, 0.8),
                force_explode=(0.3, 0.6),
                ensure_inside_axes=True,
                max_move=(20, 20),
            )
        except Exception as exc:
            warnings.warn(
                f"make_aggregated_scatter: adjust_text failed ({exc}); "
                "dropping labels rather than overlapping.",
                stacklevel=2,
            )
            for t in texts:
                t.remove()

    fig.tight_layout()
    return fig
