# ============================================================
# app.py  —  Flask backend for RegSeqDB
# ============================================================
from flask import Flask, request, render_template
from regseqdb import RegSeqDB
from plotting import (
    plot_tf_affinity_vs_expression,
    plot_rnap_energy_vs_expression,
    plot_tf_affinity_vs_rnap_energy,
    fig_to_json,
)
import pandas as pd

app = Flask(__name__)
app.config['APPLICATION_ROOT'] = '/students_26/Team12/project/app'

# ── DB credentials ────────────────────────────────────────────
DB_CREDENTIALS = dict(
    host     = "bioed-new.bu.edu",
    port     = 4253,
    database = "Team12",
    username = "mgdouq",
    password = "ha3563douq",
)

db = RegSeqDB()
db.connect(**DB_CREDENTIALS)


def get_db():
    """Return DB connection, reconnecting if it has gone away."""
    global db
    try:
        db.cursor.execute("SELECT 1")
    except Exception:
        db = RegSeqDB()
        db.connect(**DB_CREDENTIALS)
    return db


# ══════════════════════════════════════════════════════════════
# Static page routes
# ══════════════════════════════════════════════════════════════

@app.route('/')
def home():
    return render_template('defaultpage.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/help')
def help_page():
    return render_template('help.html')

@app.route('/comparison')
def comparison():
    return render_template('comparison.html')

@app.route('/singlesearch')
def singlesearch():
    return render_template('singlesearch.html')


# ══════════════════════════════════════════════════════════════
# Single search route
# ══════════════════════════════════════════════════════════════

@app.route('/search', methods=['GET'])
def search():

    promoter  = request.args.get("promoter_name", "").strip()
    tf        = request.args.get("TF_name", "").strip()
    condition = request.args.get("Condition", "").strip()

    if not promoter or not tf or not condition:
        return render_template("singlesearch.html")

    try:
        data = get_db().get_promoter_expr_and_binding(
            promoter=promoter,
            condition=condition,
            tf=tf,
            include_rnap=True,
        )
        results  = data["results"]
        colnames = data["colnames"]
        rowcount = data["rowcount"]
    except ValueError as e:
        return render_template("singlesearch.html", error=str(e))

    graph_tf_expr   = None
    graph_rnap_expr = None
    graph_tf_rnap   = None

    if rowcount > 0:
        graph_tf_expr   = fig_to_json(plot_tf_affinity_vs_expression(results, colnames))
        graph_rnap_expr = fig_to_json(plot_rnap_energy_vs_expression(results, colnames))
        graph_tf_rnap   = fig_to_json(plot_tf_affinity_vs_rnap_energy(results, colnames))

    col_idx = {name: i for i, name in enumerate(colnames)}
    table_rows = []
    for row in results:
        num_dna    = row[col_idx["num_DNA"]]
        num_rna    = row[col_idx["num_RNA"]]
        energy     = row[col_idx.get("energy",   -1)]
        affinity   = row[col_idx.get("affinity", -1)]
        expression = (num_rna / num_dna) if num_dna > 0 else 0
        table_rows.append({
            "sID":        row[col_idx["sID"]],
            "seq":        row[col_idx["seq"]] if "seq" in col_idx else "—",
            "num_DNA":    num_dna,
            "num_RNA":    num_rna,
            "expression": expression,
            "affinity":   affinity if affinity is not None else 0,
            "energy":     energy   if energy   is not None else 0,
        })

    return render_template(
        "singlesearch.html",
        promoter_name   = promoter,
        tf_name         = tf,
        condition       = condition,
        rowcount        = rowcount,
        graph_tf_expr   = graph_tf_expr,
        graph_rnap_expr = graph_rnap_expr,
        graph_tf_rnap   = graph_tf_rnap,
        table_rows      = table_rows,
        sigma_factor    = "σ70",
    )


# ══════════════════════════════════════════════════════════════
# Comparison route
# ══════════════════════════════════════════════════════════════

@app.route('/compare', methods=['GET'])
def compare():
    promoter   = request.args.get("promoter_name", "").strip()
    tf         = request.args.get("TF_name", "").strip()
    condition1 = request.args.get("Condition1", "").strip()
    condition2 = request.args.get("Condition2", "").strip()

    if not promoter or not tf or not condition1 or not condition2:
        return render_template("comparison.html")

    try:
        data1 = get_db().get_promoter_expr_and_binding(
            promoter=promoter, condition=condition1, tf=tf, include_rnap=True,
        )
        data2 = get_db().get_promoter_expr_and_binding(
            promoter=promoter, condition=condition2, tf=tf, include_rnap=True,
        )
    except ValueError as e:
        return render_template("comparison.html", error=str(e))

    results1, colnames1, rowcount1 = data1["results"], data1["colnames"], data1["rowcount"]
    results2, colnames2, rowcount2 = data2["results"], data2["colnames"], data2["rowcount"]

    graph_tf_expr   = None
    graph_rnap_expr = None
    graph_tf_rnap   = None

    if rowcount1 > 0 or rowcount2 > 0:

        def label_results(results, colnames, label):
            df = pd.DataFrame(results, columns=colnames)
            df["cond"] = label
            return df

        df_combined = pd.concat([
            label_results(results1, colnames1, condition1),
            label_results(results2, colnames2, condition2),
        ], ignore_index=True)

        combined_results = df_combined.values.tolist()
        combined_cols    = list(df_combined.columns)

        graph_tf_expr   = fig_to_json(
            plot_tf_affinity_vs_expression(
                combined_results, combined_cols,
                title=f"TF Affinity vs Expression — {condition1} vs {condition2}",
            )
        )
        graph_rnap_expr = fig_to_json(
            plot_rnap_energy_vs_expression(
                combined_results, combined_cols,
                title=f"RNAP Binding Energy vs Expression — {condition1} vs {condition2}",
            )
        )
        graph_tf_rnap   = fig_to_json(
            plot_tf_affinity_vs_rnap_energy(
                combined_results, combined_cols,
                title=f"TF Affinity vs RNAP Energy — {condition1} vs {condition2}",
            )
        )

    return render_template(
        "comparison.html",
        promoter_name   = promoter,
        tf_name         = tf,
        condition1      = condition1,
        condition2      = condition2,
        rowcount1       = rowcount1,
        rowcount2       = rowcount2,
        graph_tf_expr   = graph_tf_expr,
        graph_rnap_expr = graph_rnap_expr,
        graph_tf_rnap   = graph_tf_rnap,
    )


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5001, debug=True)