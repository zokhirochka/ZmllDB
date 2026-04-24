import json
import math
import decimal
from flask import Flask, render_template, request
from regseqdb import RegSeqDB

app = Flask(__name__)

# ── Decimal JSON Encoder ──────────────────────────────────────────
class DecimalEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, decimal.Decimal):
            return float(obj)
        return super().default(obj)

# ── Database credentials ──────────────────────────────────────────
DB_CREDENTIALS = dict(
    host     = "bioed-new.bu.edu",
    port     = 4253,
    database = "Team12",
    username = "mgdouq",
    password = "ha3563douq",
)

# ── Connect once at startup ───────────────────────────────────────
db = RegSeqDB()
try:
    db.connect(**DB_CREDENTIALS)
    print("✓ Connected to database.")
except Exception as e:
    print(f"✗ Database connection failed: {e}")


# ── Routes ────────────────────────────────────────────────────────

@app.route('/')
def home():
    return render_template('defaultpage.html')


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/help')
def help_page():
    return render_template('help.html')


@app.route('/singlesearch')
def singlesearch():
    return render_template('singlesearch.html')


@app.route('/comparison')
def comparison():
    return render_template('comparison.html')


# ── Debug ─────────────────────────────────────────────────────────
@app.route('/debug-methods')
def debug_methods():
    from regseqdb import RegSeqDB
    db_temp = RegSeqDB()
    methods = [m for m in dir(db_temp) if not m.startswith('_')]
    return f"<pre>{methods}</pre>"


# ── Single Search ─────────────────────────────────────────────────
@app.route('/search')
def search():
    promoter  = request.args.get('promoter_name', '').strip()
    tf        = request.args.get('TF_name', '').strip()
    condition = request.args.get('Condition', 'glucose').strip()

    if not promoter or not tf:
        return render_template('singlesearch.html')

    try:
        data = db.get_promoter_expr_and_binding(promoter, condition, tf)
    except ValueError as e:
        return render_template('singlesearch.html',
                               error=str(e),
                               promoter_name=promoter,
                               tf_name=tf,
                               condition=condition)

    results  = data['results']   # (sID, num_DNA, num_RNA, energy, affinity)
    rowcount = data['rowcount']

    plot_data  = []
    table_rows = []

    for row in results:
        sID, num_DNA, num_RNA, energy, affinity = row

        if num_DNA is None or num_RNA is None or num_DNA <= 0:
            continue

        plot_data.append({
            "sID":      sID,
            "dna":      float(num_DNA),
            "rna":      float(num_RNA),
            "energy":   float(energy)   if energy   is not None else None,
            "affinity": float(affinity) if affinity is not None else None,
        })

        expr = float(num_RNA) / float(num_DNA)
        table_rows.append({
            "sID":        sID,
            "seq":        "",
            "num_DNA":    num_DNA,
            "num_RNA":    num_RNA,
            "expression": expr,
            "affinity":   float(affinity) if affinity is not None else 0,
            "energy":     float(energy)   if energy   is not None else 0,
        })

    has_data   = len(plot_data) > 0
    graph_json = json.dumps(plot_data, cls=DecimalEncoder) if has_data else None

    return render_template('singlesearch.html',
        promoter_name  = promoter,
        tf_name        = tf,
        condition      = condition,
        rowcount       = rowcount,
        graph_json     = graph_json,
        graph_tf_expr  = has_data,
        graph_rnap_expr= has_data,
        graph_tf_rnap  = has_data,
        table_rows     = table_rows,
    )


# ── Comparison ────────────────────────────────────────────────────
@app.route('/compare')
def compare():
    promoter   = request.args.get('promoter_name', '').strip()
    tf         = request.args.get('TF_name', '').strip()
    condition1 = request.args.get('Condition1', 'glucose').strip()
    condition2 = request.args.get('Condition2', 'xylose').strip()

    if not promoter or not tf:
        return render_template('comparison.html')

    try:
        inputs = [promoter, condition1, tf, promoter, condition2]
        query = """
            SELECT sID, e1.num_DNA as Cond1_DNA, e1.num_RNA as Cond1_RNA,
                   e2.num_DNA as Cond2_DNA, e2.num_RNA as Cond2_RNA, energy, affinity
            FROM (SELECT sID, SUM(num_DNA) as num_DNA, SUM(num_RNA) as num_RNA, energy, affinity
                  FROM Promoters
                  JOIN PromoterSequences AS ps USING(pID)
                  JOIN BarcodeCounts USING(sID)
                  JOIN Experiments USING (eID)
                  JOIN BindingSitesRNAP using(sID)
                  JOIN BindingSitesTF USING(sID)
                  JOIN TranscriptionFactors USING(tID)
                  WHERE pID = (SELECT pID FROM Promoters WHERE pro_name = %s)
                        AND eID = (SELECT eID FROM Experiments WHERE cond = %s LIMIT 1)
                        AND tID = (SELECT tID FROM TranscriptionFactors WHERE tf_name = %s)
                  GROUP BY ps.seq) as e1
            JOIN (SELECT sID, SUM(num_DNA) as num_DNA, SUM(num_RNA) as num_RNA
                  FROM Promoters
                  JOIN PromoterSequences AS ps USING(pID)
                  JOIN BarcodeCounts USING(sID)
                  JOIN Experiments USING (eID)
                  WHERE pID = (SELECT pID FROM Promoters WHERE pro_name = %s)
                        AND eID = (SELECT eID FROM Experiments WHERE cond = %s LIMIT 1)
                  GROUP BY ps.seq) as e2
            USING (sID);
        """
        import utils
        results_dict = utils.exec_query(cursor=db.cursor, query=query, inputs=inputs)
        results  = results_dict['results']
        rowcount = results_dict['rowcount']
    except Exception as e:
        return f"<pre>QUERY ERROR:\n{e}</pre>", 500

    try:
        plot_data = []
        for row in results:
            sID, dna1, rna1, dna2, rna2, energy, affinity = row
            dna1 = float(dna1) if dna1 is not None else 0
            rna1 = float(rna1) if rna1 is not None else 0
            dna2 = float(dna2) if dna2 is not None else 0
            rna2 = float(rna2) if rna2 is not None else 0
            if dna1 > 0 and rna1 > 0 and dna2 > 0 and rna2 > 0:
                plot_data.append({
                    "sID":      sID,
                    "dna1":     dna1,
                    "rna1":     rna1,
                    "dna2":     dna2,
                    "rna2":     rna2,
                    "energy":   float(energy)   if energy   is not None else None,
                    "affinity": float(affinity) if affinity is not None else None,
                })
        has_data = len(plot_data) > 0
        graph_payload = {"rows": plot_data, "cond1": condition1, "cond2": condition2}
        graph_json = json.dumps(graph_payload, cls=DecimalEncoder) if has_data else None
    except Exception as e:
        return f"<pre>DATA PROCESSING ERROR:\n{e}</pre>", 500

    return render_template('comparison.html',
        promoter_name=promoter, tf_name=tf,
        condition1=condition1, condition2=condition2,
        rowcount1=rowcount, rowcount2=rowcount,
        graph_tf_expr=has_data, graph_json=graph_json,
    )


# ── Run ───────────────────────────────────────────────────────────
if __name__ == '__main__':
    app.run(debug=True)
