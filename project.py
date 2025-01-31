from flask import Flask, render_template, request
from flask_sqlalchemy import SQLAlchemy
import os

app = Flask(__name__)
# below give link to database on cloud
app.config['SQLALCHEMY_DATABASE_URI'] = f"sqlite:///{os.path.expanduser('~/Downloads/data.db')}"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class SNP(db.Model):
    # id = db.Column(db.Integer)
    SNPS = db.Column(db.String(100), unique=True, primary_key=True)
    MAPPED_GENE = db.Column(db.String(100), nullable=True)
    CHR_POS = db.Column(db.Integer, unique=True)

class MappedGene(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), unique=True, nullable=False)
    ontology_term = db.Column(db.String(100), nullable=True)

class PopulationData(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    population_name = db.Column(db.String(100), nullable=False)
    statistics = db.Column(db.String(100), nullable=True)

# home page
@app.route('/')
def index():
    return render_template("Home_page.html")

# SNP query page
@app.route('/query/', methods=['GET', 'POST'])
def query():
    if request.method == 'POST':
        # grab data from 3 fields
        query1 = request.form['query1']
        query2 = request.form['query2']
        query3 = request.form['query3']
        # query the SQL file, and filter to specifications
        snps = SNP.query.filter((SNP.SNPS == query1) | (SNP.MAPPED_GENE == query2) | (SNP.CHR_POS == query3)).all()
        # return the query output
        return render_template("query_result.html", snps=snps)

    # Render query form, "GET" method
    return render_template("query.html")


# ignore the code below, still in production
@app.route('/gene/<gene_name>')
def gene_details(gene_name):
    gene_info = Gene.query.filter_by(name=gene_name).first()
    return render_template('gene_details.html', gene=gene_info)

@app.route('/population', methods=['GET', 'POST'])
def population():
    if request.method == 'POST':
        population = request.form['population']
        stats = PopulationData.query.filter_by(population_name=population).all()
        return render_template('population_stats.html', stats=stats)
    return render_template('population.html')

@app.route('/visualization')
def visualization():
    return render_template('visualization.html')

@app.route('/api/snp/<rsid>', methods=['GET'])
def get_snp(rsid):
    snp = SNP.query.filter_by(rsid=rsid).first()
    if snp:
        return jsonify({"rsid": snp.rsid, "gene_name": snp.gene_name, "position": snp.position, "p_value": snp.p_value})
    return jsonify({"error": "SNP not found"}), 404


if __name__ == '__main__':
    app.run(debug=True)

# test comment
