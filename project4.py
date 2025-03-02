from flask import Flask, jsonify, render_template, request, send_file
from flask_sqlalchemy import SQLAlchemy
import pandas as pd
import os
from io import BytesIO
import base64
import plotly.express as px 

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


class NSL(db.Model):
    __tablename__ = 'Combined_Table'
    id = db.Column('rowid', db.Integer, primary_key=True)
    chrom_num = db.Column('Chromosome_Number', db.Integer, unique=False)
    win_start = db.Column('Window_Size_Lower', db.Integer)
    nsl_val = db.Column('Standardised_nSL', db.Float, nullable=True)
    nsl_sa_pop = db.Column('population', db.String(50))
    nsl_pos = db.Column('Physical_Position', db.Integer)
    


# home page
@app.route('/')
def index():
    return render_template("index.html")

# about page
@app.route('/about')
def about():
    return render_template("about.html")


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
   
    # get data from input rsid for pos slection
    gene = snp.mapped_gene
    p_value = snp.snp_p_value
    phenotype = snp.snp_phenotype
    population = snp.snp_population
    chromosome = snp.chr_id
    position = snp.gene_pos
    
    # get the query info
    if request.method == 'POST':
        query5 = request.form.get('query5', '')
        query6 = request.form.get('query6', '')
        query7 = request.form.get('query7', '')
         
        # display the population info
        pop_info_disp = population_info.get(query5, "Please select a population")
         
        # if Tajimas or nSL stat is selected...
        if query6 in ["Tajima's D", "nSL"]:

            if query6 == "Tajima's D":
                StatModel = Tajima
                columns = {
                    "chrom": "chrom",
                    "bin_start": "bin_start",
                    "stat_column": "tajD",
                    "population": "sa_pop",
                    "plot_title": "Tajima's D"}
                
                stat_title = "Tajima's D"

                # Calculate window boundaries
                window = (position // 10000) * 10000
                lower = window - (int(query7) * 1000)
                upper = window + (int(query7) * 1000) 
            
                # Query the DB
                if query5 == "All":
                    filt = db.session.query(StatModel).filter(
                        getattr(StatModel, columns['chrom']) == chromosome,
                        getattr(StatModel, columns['bin_start']).between(lower, upper)
                    ).all()
                else:
                    filt = db.session.query(StatModel).filter(
                        getattr(StatModel, columns['population']) == query5,
                        getattr(StatModel, columns['chrom']) == chromosome,
                        getattr(StatModel, columns['bin_start']).between(lower, upper)
                    ).all()


                if filt:
                    df = pd.DataFrame([row.__dict__ for row in filt])
                    df.drop(columns=['_sa_instance_state'], inplace=True)

                    # Compute region stats
                    region_mean = round(df[columns['stat_column']].mean(), 4)
                    region_std = round(df[columns['stat_column']].std(), 4)
                    stat_value = next(
                        (getattr(row, columns['stat_column']) 
                        for row in filt 
                        if getattr(row, columns['bin_start']) == window),
                        None
                    )

                    df[f"{columns['plot_title']} Mean (Selected Region)"] = region_mean
                    df[f"{columns['plot_title']} Std. (Selected Region)"] = region_std
                    df[columns['plot_title']] = stat_value

                    # to show the tajima d value to user
                    final_stat = df.iloc[0,7]

                    # --- PLOTLY CODE ---
                    # Build an interactive scatter plot
                    # color by population if user selected "All"
                    color_col = columns['population'] if query5 == "All" else None

                    fig = px.scatter(
                        df,
                        x=columns['bin_start'],
                        y=columns['stat_column'],
                        color=color_col,  # color by population only if "All"
                        title=f"{columns['plot_title']} for Chromosome Window {window} ± {query7} kb",
                        hover_data=[columns['bin_start'], columns['stat_column'], columns['population']]
                    )
                    fig.update_layout(
                        xaxis_title=f"Chromosome {chromosome} Region (bp)",
                        yaxis_title=columns['plot_title']
                    )

                    # Add horizontal lines (threshold, region mean)
                    # Plotly's add_hline is in plotly.graph_objects.Figure
                    # but we can do it via fig.add_hline():
                    fig.add_hline(y=-2, line_color="red", annotation_text="Threshold (-2)", annotation_position="bottom left")
                    fig.add_hline(y=region_mean, line_color="green", annotation_text="Region Mean", annotation_position="bottom left")

                    # Convert to HTML for rendering in the template
                    plot_html = fig.to_html(full_html=False)

                    # Handle the "Download TSV" request
                    if 'download' in request.form:
                        tsv_buffer = BytesIO()
                        df.to_csv(tsv_buffer, index=False, sep="\t")
                        tsv_buffer.seek(0)
                        return send_file(
                            tsv_buffer,
                            mimetype="text/tab-separated-values",
                            as_attachment=True,
                            download_name=f"{columns['plot_title']}_{rs_value}.tsv"
                        )

                    # Return the Plotly HTML in the template
                    return render_template(
                        'visualisation2.html',  # <-- A separate template or rename existing
                        rs_value=rs_value,
                        plot_html=plot_html,
                        snp=snp,
                        pop_info_disp=pop_info_disp,
                        filt=filt,
                        region_mean=region_mean,
                        region_std=region_std,
                        stat_value=stat_value,
                        stat_title=stat_title,
                        final_stat=final_stat
                    )



            elif query6 == "nSL":
                StatModel = NSL
                columns = {
                    "chrom": "chrom_num",
                    "bin_start": "win_start",
                    "stat_column": "nsl_val",
                    "population": "nsl_sa_pop",
                    "plot_title": "nSL",
                    "nsl_pos" : "nsl_pos"}
                
                stat_title = "nSL"

                # Calculate window boundaries
                lower = position - (int(query7) * 1000) 
                upper = position + (int(query7) * 1000)
                window = position

                # Query the DB
                if query5 == "All":
                    filt = db.session.query(StatModel).filter(
                        getattr(StatModel, columns['chrom']) == chromosome,
                        NSL.nsl_pos.between(lower, upper)
                    ).all()
                else:
                    filt = db.session.query(StatModel).filter(
                        getattr(StatModel, columns['population']) == query5,
                        getattr(StatModel, columns['chrom']) == chromosome,
                        NSL.nsl_pos.between(lower, upper)
                    ).all()

                if filt:
                    df = pd.DataFrame([row.__dict__ for row in filt])
                    df.drop(columns=['_sa_instance_state'], inplace=True)

                    # Compute region stats
                    df["search_position"] = position
                    region_mean = round(df[columns['stat_column']].mean(), 4)
                    region_std = round(df[columns['stat_column']].std(), 4)
                    stat_value = next(
                        (getattr(row, columns['stat_column']) 
                        for row in filt 
                        if getattr(row, columns['bin_start']) == window),
                        None
                    )

                    df[f"{columns['plot_title']} Mean (Selected Region)"] = region_mean
                    df[f"{columns['plot_title']} Std. (Selected Region)"] = region_std
                    df[columns['plot_title']] = stat_value
                    
                    # present nSL value to user
                    # final_stat = df[df['nsl_pos'] == position]['nsl_val']
                    matched_rows = df[df['search_position'] == df['nsl_pos']]

                    print(df)
                    print(matched_rows.head())
                    # Get the 'nsl_val' from the first matched row
                    final_stat = matched_rows.iloc[0]['nsl_val'] if not matched_rows.empty else None


                    # PLOTLY CODE
                    # Build an interactive scatter plot
                    # color by population if user selected "All"
                    color_col = columns['population'] if query5 == "All" else None

                    fig = px.scatter(
                        df,
                        x=columns['nsl_pos'],
                        y=columns['stat_column'],
                        color=color_col,  # color by population only if "All"
                        title=f"{columns['plot_title']} for Chromosome Positon {position} ± {query7} kb",
                        hover_data=[columns['bin_start'], columns['stat_column'], columns['population']]
                    )
                    fig.update_layout(
                        xaxis_title=f"Chromosome {chromosome} Region (bp)",
                        yaxis_title=columns['plot_title'],
                        xaxis=dict(
                        range=[lower, upper],
                        showgrid=True,
                        zeroline=True,
                        )
                    )

                    # Add horizontal lines (threshold, region mean)
                    # Plotly's add_hline is in plotly.graph_objects.Figure
                    # but we can do it via fig.add_hline():
                    fig.add_hline(y=-2, line_color="red", annotation_text="Threshold (-2)", annotation_position="bottom left")
                    fig.add_hline(y=region_mean, line_color="green", annotation_text="Region Mean", annotation_position="bottom left")
                    fig.add_vline(x=position, line=dict(color="black", width=2, dash="dash"))

                    # Convert to HTML for rendering in the template
                    plot_html = fig.to_html(full_html=False)

                    # Handle the "Download TSV" request
                    if 'download' in request.form:
                        tsv_buffer = BytesIO()
                        df.to_csv(tsv_buffer, index=False, sep="\t")
                        tsv_buffer.seek(0)
                        return send_file(
                            tsv_buffer,
                            mimetype="text/tab-separated-values",
                            as_attachment=True,
                            download_name=f"{columns['plot_title']}_{rs_value}.tsv"
                        )

                    # Return the Plotly HTML in the template
                    return render_template(
                        'visualisation2.html',  
                        rs_value=rs_value,
                        plot_html=plot_html,
                        snp=snp,
                        pop_info_disp=pop_info_disp,
                        filt=filt,
                        region_mean=region_mean,
                        region_std=region_std,
                        stat_value=stat_value,
                        stat_title=stat_title,
                        final_stat=final_stat
                    )

    # Default GET request or no data
    return render_template('visualisation2.html', rs_value=rs_value, snp=snp)


@app.route('/query/visualisation/<rs_value>/download', methods=['GET', 'POST'])
def download_stats(rs_value):
    pop = request.args.get('pop', 'All')
    window_size = int(request.args.get('window', 10))  # Default 10kb window

    snp = db.session.query(SNP).filter_by(rs_value=rs_value).first()
    if not snp:
        return "SNP not found", 404

    chromosome = snp.chr_id
    position = snp.gene_pos

    window = (position // 10000) * 10000
    lower = window - (window_size * 1000)
    upper = window + (window_size * 1000)

    query = db.session.query(Tajima).filter(
        Tajima.chrom == chromosome,
        Tajima.bin_start.between(lower, upper)
    )
    
    if pop != "All":
        query = query.filter(Tajima.sa_pop == pop)

    filt = query.all()

    if not filt:
        return "No Tajima's D data found for this SNP.", 404

    df = pd.DataFrame([row.__dict__ for row in filt])
    df.drop(columns=['_sa_instance_state'], inplace=True)

    #average, std and SNPs Tajimas 
    region_mean = round(df['tajD'].mean(), 4)
    region_std = round(df['tajD'].std(),4)
    tajD_value= next((row.tajD for row in filt if row.bin_start == window))
    #now adding these calulations to named colums in the df
    df["Tajima's D Mean (Selected Region)"] = region_mean
    df["Tajima's D Std. (Selected Region)"] = region_std   
    df["Tajima's D"] = tajD_value

    tsv_buffer = BytesIO()
    df.to_csv(tsv_buffer, index=False, sep="\t")
    tsv_buffer.seek(0)

    return send_file(
        tsv_buffer,
        mimetype="text/tab-separated-values",
        as_attachment=True,
        download_name=f"tajimasD_{rs_value}.tsv"
    )



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
