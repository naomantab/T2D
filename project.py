from flask import Flask, jsonify, render_template, request
from flask_sqlalchemy import SQLAlchemy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
from io import BytesIO
import base64


app = Flask(__name__)
# below give link to database on cloud
app.config['SQLALCHEMY_DATABASE_URI'] = f"sqlite:///{os.path.expanduser('~/Downloads/t2d_database.db')}"
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

class Tajima(db.Model):
    __tablename__ = 'TajimasD_results_ALL_POPULATIONS' 
    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column('CHROM', db.Integer, unique=False)
    bin_start = db.Column('BIN_START', db.Integer, unique= False)
    tajD = db.Column('TajimaD', db.Float, nullable=True)
    sa_pop = db.Column('Population', db.String(50))

#class iHS(db.Model):
    #__tablename__ = ''


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


@app.route('/query/visualisation/<rs_value>/', methods=['GET', 'POST'])
def visualisation(rs_value):
    snp = db.session.query(SNP).filter_by(rs_value=rs_value).first()
   
   #get data from input rsid for pos slection
    gene = snp.mapped_gene
    p_value = snp.snp_p_value
    phenotype = snp.snp_phenotype
    population = snp.snp_population
    chromosome= snp.chr_id
    position= snp.gene_pos
    #image= 

    

    

    if request.method == 'POST':
         query5 = request.form.get('query5', '')
         query6 = request.form.get('query6', '')
         query7 = request.form.get('query7', '')
         
         pop_info_disp = population_info.get(query5, "Please select a population")
         
         if query6 == "Tajima's D":
             window= (position // 10000) * 10000
             #make kb
             lower= window - (int(query7) * 1000)
             upper= window + (int(query7) * 1000)

             if query5 == "All":
                filt = db.session.query(Tajima).filter(
                    Tajima.chrom == chromosome,
                    Tajima.bin_start.between(lower, upper)
                    ).all()
             else:
                filt = db.session.query(Tajima).filter(
                    Tajima.sa_pop == query5,
                    Tajima.chrom == chromosome,
                    Tajima.bin_start.between(lower, upper)
                    ).all()
                   
             if filt:
                df = pd.DataFrame([row.__dict__ for row in filt])  
                df.drop(columns=['_sa_instance_state'], inplace=True)
                   
    
                #print(df)
            
                #average and st of region of inrest
                region_mean = round(df['tajD'].mean(), 4)
                region_std = round(df['tajD'].std(),4)

                #clear previous plot just in case
                plt.clf()
                plt.figure(figsize=(10,6))
                tajD_value= next((row.tajD for row in filt if row.bin_start == window))
                #plot all or plot 1 population
                if query5 == "All":
                    populations = df['sa_pop'].unique()
                    colors = plt.cm.get_cmap('tab10', len(populations))

                    for idx, population in enumerate(populations):
                        population_data = df[df['sa_pop'] == population]
                        plt.scatter(population_data['bin_start'], population_data['tajD'], alpha=0.6, 
                                    label=population, color=colors(idx))

                else:
                    plt.scatter(df['bin_start'], df['tajD'], alpha=0.6, label=query5, color='blue')


         
                #plot figure
                
                plt.axhline(y=-2, color = 'red', linestyle= '-')
                plt.xlabel(f"Chromsome {chromosome} Region (bp)")
                plt.ylabel("Tajima's D")
                plt.legend(loc='upper right')
                plt.title(f"Tajima's D for Chromosome Position {window} ± {query7} kb")
                

                #save plt 
                buf= BytesIO()
                plt.savefig(buf, format= "png")
                plt.close()

                data= base64.b64encode(buf.getbuffer()).decode("ascii")
                return render_template('visualisation.html', rs_value=rs_value, image_data= data, snp=snp, pop_info_disp=pop_info_disp, filt=filt, region_mean=region_mean, region_std=region_std, tajD_value=tajD_value)

                
        
        #if query5 and query6 == "nSL":
        ### need to add data to db to make filter query

       #nSLPlot()
   
    return render_template('visualisation.html', rs_value=rs_value, snp=snp)


# Info for each population
global population_info
population_info = {
        "Bengali in Bangladesh": """🌎 Bengali in Bangladesh. 
                                This population is sampled directly from Bangladesh, representing the genetic diversity of ethnic Bengalis living in the region. 
                                Data comes from Phase 3 of the 1000 Genomes Project, with  a total of 87 samples. 
                                These participants recruited from rural and urban areas across Bangladesh.""",
        "Gujarati Indians in Houston, TX": """🌎 Gujarati Indians in Houston, Texas. 
                                        This is a diaspora population, with 107 participants of Gujarati Indian ancestry living in Houston, Texas. 
                                        This group is part of Phase 3 of the 1000 Genomes Project. 
                                        Their genetic data is especially intruiging because of potential environmental effects linked to migration, diet changes, and lifestyle shifts in the US.""",
        "Indian Telugu in the UK": """🌎 Indian Telugu in the UK. 
                                 This group includes 103 individuals of Telugu-speaking Indian ancestry residing in the United Kingdom, many of whom migrated in the 20th century. 
                                 Data comes from Phase 3 of the 1000 Genomes Project.""",
        "Punjabi in Lahore, Pakistan": """🌎 Punjabi in Lahore, Pakistan. 
                                    This population consists of 97 ethnic Punjabis living in Lahore, Pakistan. 
                                    Sampled locally fir Phase 3 of the 1000 Genomes Project, these participants provide a local reference population for understanding genetic variation and T2D risk within South Asia.""",
        "All": "Please select a population to see the information."
}


if __name__ == '__main__':
    app.run(debug=True)

# test comment
