"""
plotting.py
-----------
Visualization module for the RegSeqDB website project.
Based on Dr. James Galagan's lab data (BU Microbiology).
"""

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
PALETTE       = ["#00bfff"]
BG_COLOR      = "#0f1117"
PAPER_COLOR   = "#1a1d27"
FONT_COLOR    = "#e8eaf0"
GRID_COLOR    = "#2e3147"

BASE_LAYOUT = dict(
    paper_bgcolor=PAPER_COLOR,
    plot_bgcolor=BG_COLOR,
    font=dict(color=FONT_COLOR, family="monospace", size=13),
    margin=dict(l=60, r=30, t=60, b=60),
    hoverlabel=dict(bgcolor="#2e3147", font_size=12),
    xaxis=dict(gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR),
    yaxis=dict(gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _to_df(results, colnames):
    return pd.DataFrame(results, columns=colnames)


def _compute_expression(df):
    df = df[df["num_DNA"] > 0].copy()
    df["expression"] = df["num_RNA"] / df["num_DNA"]
    df = df[df["expression"] > 0]
    return df


def _apply_log_y(fig, enabled):
    if enabled:
        fig.update_yaxes(type="log", gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR)
    return fig


# ---------------------------------------------------------------------------
# Plot 1 — TF Affinity vs. Expression
# ---------------------------------------------------------------------------
def plot_tf_affinity_vs_expression(results, colnames, title=None, log_scale=True):
    df = _to_df(results, colnames)
    df = _compute_expression(df)

    color_col = "tf_name" if "tf_name" in df.columns else None

    fig = px.scatter(
        df,
        x="affinity",
        y="expression",
        color=color_col,
        hover_data=["sID"] + (["tf_name"] if color_col else []),
        labels={
            "affinity":   "TF Binding Affinity",
            "expression": "Expression (RNA / DNA)",
            "tf_name":    "Transcription Factor",
        },
        color_discrete_sequence=PALETTE,
        opacity=0.75,
        trendline="ols",
        trendline_options=dict(log_y=log_scale),
    )

    fig.update_traces(marker=dict(size=7, line=dict(width=0.5, color="rgba(255,255,255,0.19)")))
    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=title or "TF Binding Affinity vs. Promoter Expression",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot 2 — RNAP Binding Energy vs. Expression
# ---------------------------------------------------------------------------
def plot_rnap_energy_vs_expression(results, colnames, title=None, log_scale=True):
    df = _to_df(results, colnames)
    df = _compute_expression(df)

    color_col = "sigma" if "sigma" in df.columns else None

    fig = px.scatter(
        df,
        x="energy",
        y="expression",
        color=color_col,
        hover_data=["sID"],
        labels={
            "energy":     "RNAP Binding Energy",
            "expression": "Expression (RNA / DNA)",
            "sigma":      "Sigma Factor",
        },
        color_discrete_sequence=PALETTE,
        opacity=0.75,
        trendline="ols",
        trendline_options=dict(log_y=log_scale),
    )

    fig.update_traces(marker=dict(size=7, symbol="diamond", line=dict(width=0.5, color="rgba(255,255,255,0.19)")))
    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=title or "RNAP Binding Energy vs. Promoter Expression",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot 3 — Expression by Condition
# ---------------------------------------------------------------------------
def plot_expression_by_condition(results, colnames, title=None, log_scale=True):
    df = _to_df(results, colnames)
    df = _compute_expression(df)

    fig = px.box(
        df,
        x="cond",
        y="expression",
        color="cond",
        points="all",
        hover_data=["sID"],
        labels={
            "cond":       "Experimental Condition",
            "expression": "Expression (RNA / DNA)",
        },
        color_discrete_sequence=PALETTE,
    )

    fig.update_traces(marker=dict(size=4, opacity=0.6))
    fig.update_layout(
        **BASE_LAYOUT,
        showlegend=False,
        title=dict(
            text=title or "Promoter Expression Across Conditions",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot 4 — TF Affinity vs. RNAP Energy
# ---------------------------------------------------------------------------
def plot_tf_affinity_vs_rnap_energy(results, colnames, title=None):
    df = _to_df(results, colnames)

    color_col  = "tf_name" if "tf_name" in df.columns else None
    symbol_col = "sigma"   if "sigma"   in df.columns else None

    fig = px.scatter(
        df,
        x="affinity",
        y="energy",
        color=color_col,
        symbol=symbol_col,
        hover_data=["sID"],
        labels={
            "affinity": "TF Binding Affinity",
            "energy":   "RNAP Binding Energy",
            "tf_name":  "Transcription Factor",
            "sigma":    "Sigma Factor",
        },
        color_discrete_sequence=PALETTE,
        opacity=0.8,
    )

    fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color="rgba(255,255,255,0.19)")))
    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=title or "TF Affinity vs. RNAP Binding Energy",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    return fig


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------
def fig_to_json(fig):
    import json
    import plotly
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
