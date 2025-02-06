from flask import Flask, jsonify, render_template, request
from flask_sqlalchemy import SQLAlchemy
import os

app = Flask(__name__)
# below give link to database on cloud
app.config['SQLALCHEMY_DATABASE_URI'] = f"sqlite:///{os.path.abspath('t2d_database.db')}"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class SNP(db.Model):
    # id = db.Column(db.Integer)
    rs_value = db.Column('SNPS',db.String(100), unique=True, primary_key=True)
    gene_pos = db.Column('CHR_POS',db.Integer, unique=True)
    mapped_gene = db.Column('MAPPED_GENE',db.String(100), nullable=True)
    snp_p_value = db.Column('P-VALUE',db.Float, nullable=True)
    snp_phenotype = db.Column('DISEASE/TRAIT',db.Float, nullable=True)
    snp_population = db.Column('INITIAL SAMPLE SIZE',db.Float, nullable=True)

class Ontology(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    ont_mapped_gene = db.Column(db.String(100), unique=True, nullable=False)
    ontology_term = db.Column(db.String(100), nullable=True)

class PopulationData(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    population_name = db.Column(db.String(100), nullable=False)
    statistics = db.Column(db.String(100), nullable=True)

# home page
@app.route('/')
def index():
    return render_template("index.html")

# SNP query page
@app.route('/query/', methods=['GET', 'POST'])
def query():
    snps = None  # Default value when the page is first loaded
    if request.method == 'POST':
        query1 = request.form.get('query1', '')
        query2 = request.form.get('query2', '')
        query3 = request.form.get('query3', '')
        # Search query only if at least one input is filled
        snps = SNP.query.filter((SNP.rs_value == query1) | (SNP.gene_pos == query2) | (SNP.mapped_gene == query3)).all()
    return render_template("query.html", snps=snps)


# ignore the code below, still in production
@app.route('/query/ontology/<mapped_gene>/')
def ontology(mapped_gene):
    # ontology_info = ontology.query.filter_by(name=gene_name).first()
    return render_template('ontology.html', ont_mapped_gene = mapped_gene)






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
