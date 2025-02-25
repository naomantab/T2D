from flask import Flask, jsonify, render_template, request
from flask_sqlalchemy import SQLAlchemy
import os

app = Flask(__name__)
# below give link to database on cloud
app.config['SQLALCHEMY_DATABASE_URI'] = f"sqlite:///{os.path.expanduser('~/Downloads/TAJD_STATS.db')}"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)


class SNP(db.Model):
    __tablename__ = 'SNP_ALL'
    rs_value = db.Column('SNPS', db.String(100), unique=True, primary_key=True)
    gene_pos = db.Column('CHR_POS', db.Integer, unique=True)
    mapped_gene = db.Column('MAPPED_GENE', db.String(100), nullable=True)
    snp_p_value = db.Column('P-VALUE', db.Float, nullable=True)
    snp_phenotype = db.Column('MAPPED_TRAIT', db.Float, nullable=True)
    snp_population = db.Column('General Ancestry', db.Float, nullable=True)
    chr_id = db.Column('CHR_ID', db.Integer, nullable=True)


class Ontology(db.Model):
    __tablename__ = 'ontology_term'
    go_id = db.Column('GO_ID', db.String(50), nullable=True, primary_key=True)
    gene_name = db.Column('Gene_Name', db.String(100), nullable=True)
    qualifier = db.Column('Qualifier', db.String(50), nullable=True)
    gene_function = db.Column('Gene_Function', db.String(255), nullable=True)
    evidence_code = db.Column('Evidence_Code', db.String(50), nullable=True)
    aspect = db.Column('Aspect', db.String(5), nullable=True)


class PopulationData(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    population_name = db.Column(db.String(100), nullable=False)
    statistics = db.Column(db.String(100), nullable=True)

class SNP_sumstat(db.Model):
    __tablename__ = 'SNP_STATS'
    mapped_gene_stats = db.Column('MAPPED_GENE', db.String(100), nullable=True, primary_key=True)
    risk_mean = db.Column('MEAN_RISK_ALLELE_FREQUENCY', db.Float, nullable=True)
    risk_std = db.Column('STD_RISK_ALLELE_FREQUENCY', db.Float, nullable=True)

class plot(db.Model):
    __tablename__ = 'TAJD_BEB'
    id = db.Column(db.Integer, primary_key=True) 
    chrom = db.Column('CHROM', db.Integer, unique=False)
    bin_start = db.Column('BIN_START', db.Integer, unique= False)
    tajD = db.Column('TajimaD', db.Integer, nullable=True)
    sa_pop = db.Column('Population', db.String(50))

# home page


@app.route('/')
def index():
    return render_template("index.html")

# SNP query page


@app.route('/query/', methods=['GET', 'POST'])
def query():
    snps = None  # Default value when the page is first loaded
    if request.method == 'POST':
        # Grab user input
        query1 = request.form.get('query1', '')
        query2 = request.form.get('query2', '')
        query3 = request.form.get('query3', '')
        query4 = request.form.get('query4', '')
            
        # Join the 2 tables and query using the input
        if query1 or query2 or query3 or query4:
            snps = db.session.query(SNP, SNP_sumstat).join(SNP_sumstat, SNP.mapped_gene == SNP_sumstat.mapped_gene_stats).filter(
                (SNP.rs_value == query1) | 
                (SNP.gene_pos == query2) | 
                (SNP.snp_population == query4) |
                (SNP.mapped_gene == query3)).all()
                
    return render_template("query.html", snps=snps)


# ignore the code below, still in production
@app.route('/query/ontology/<mapped_gene>/')
def ontology(mapped_gene):
    ontology_results = Ontology.query.filter_by(gene_name=mapped_gene).all()

    # Render the ontology results page
    return render_template('ontology.html', mapped_gene=mapped_gene, ontology_results=ontology_results)

@app.route('/population', methods=['GET', 'POST'])
def population():
    if request.method == 'POST':
        population = request.form['population']
        stats = PopulationData.query.filter_by(
            population_name=population).all()
        return render_template('population_stats.html', stats=stats)
    return render_template('population.html')


@app.route('/query/visualization/<rs_value>/', methods=['GET', 'POST'])
def visualization(rs_value):
    snps = SNP.query.filter_by(rs_value=rs_value).first()  # Get SNP data for the given rs_value
    
    if snps:  # Check if the SNP record exists
        user_chromosome = snps.chr_id  # Retrieve the chromosome from SNP record
    
        if request.method == 'POST':
            query5 = request.form.get('query5', '')
            query6 = request.form.get('query6', '')

            if query5 and query6 == "tajD":
                snp_record = SNP.query.filter_by(rs_value=rs_value).first()

                if snp_record:
                    target = snp_record.gene_pos
                    region_size = 10000

                # Query the plot table for chromosome data related to the gene
                    plot_data = plot.query.filter_by(chrom=user_chromosome).all()  # Assuming chrom=3 as the relevant chromosome

                    # Convert the result to a list of BIN_START values for computation
                    bin_start_values = [entry.bin_start for entry in plot_data]

                    # Find the closest BIN_START to the SNP gene position
                    closest_wndw = min(bin_start_values, key=lambda x: abs(x - target))

                    # Filter the plot data for regions around this closest BIN_START
                    filtered_plot_data = [entry for entry in plot_data if 
                                        closest_wndw - region_size <= entry.bin_start <= closest_wndw + region_size]

                    # Now you can use `filtered_plot_data` for your visualization
                    # Pass the data to the template or generate the plot as needed
                    return render_template('visualization.html', plot_data=filtered_plot_data, closest_wndw=closest_wndw)
                
                #  plot_q = db.session.query(
                #     (plot.sa_pop == query5) & 
                #     (plot.bin_start == (snps.gene_pos )).all() 
    return render_template('visualization.html')


def TajDPlot(plot_data, closest_wndw, region_size):
    # Prepare data for the plot
    x_vals = [entry.bin_start for entry in plot_data]
    y_vals = [entry.tajD for entry in plot_data]

    # Plotting the data
    plt.scatter(x_vals, y_vals, alpha=0.6)
    plt.axhline(y=-2, color='red', linestyle='-')  # Tajima's D threshold line
    plt.xlabel(f"Chromosome 3 Region (bp)")
    plt.ylabel("Tajima's D")
    plt.title(f"Tajima's D for Chromosome Position {closest_wndw} ± {region_size}")
    
    # Show or save the plot
    plt.show()  # or plt.savefig("tajd_plot.png")

if __name__ == '__main__':
    app.run(debug=True)

# test comment
