"""
Flask App for T2D Analysis

This application serves a web interface to query and visualize
SNP data, including Tajima's D and nSL statistics.

Authors: Naoman Tabassam, Deepak Raj, Katelin Cunningham, Olayemi Bakare, Joseph Eytle
Date: March 2025
"""

# import all relevant modules 
from flask import Flask, render_template, request, send_file
from flask_sqlalchemy import SQLAlchemy
import pandas as pd
import os
from io import BytesIO
import base64
import plotly.express as px
import numpy as np
import requests
import json

# intitalise the flask app
app = Flask(__name__)

# below give link to database 
app.config['SQLALCHEMY_DATABASE_URI'] = f"sqlite:///{os.path.expanduser('~/Downloads/t2d_database.db')}"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
# initialise flask app with SQL-alchemy ORM
db = SQLAlchemy(app)

# link each table from database to classes on python

# class with SNP data with all populations
class SNP(db.Model):
    __tablename__ = 'SNP'
    rs_value = db.Column('SNP_rsID', db.String(100), unique=True, primary_key=True)
    gene_pos = db.Column('CHR_POS', db.Integer, unique=True)
    mapped_gene = db.Column('MAPPED_GENE', db.String(100), nullable=True)
    snp_p_value = db.Column('P-VALUE', db.Float, nullable=True)
    snp_phenotype = db.Column('MAPPED_TRAIT', db.String, nullable=True)
    snp_population = db.Column('GENERAL ANCESTRY', db.String, nullable=True)
    chr_id = db.Column('CHR_ID', db.Integer, nullable=True)

# class with ontology data
class Ontology(db.Model):
    __tablename__ = 'ontology_term'
    go_id = db.Column('GO_ID', db.String(50), nullable=True, primary_key=True)
    gene_name = db.Column('Gene_Name', db.String(100), nullable=True)
    qualifier = db.Column('Qualifier', db.String(50), nullable=True)
    gene_function = db.Column('Gene_Description', db.String(255), nullable=True)
    evidence_code = db.Column('Evidence_Code', db.String(50), nullable=True)
    aspect = db.Column('Aspect', db.String(5), nullable=True)

# class for gene wide summary statistics
class SNP_sumstat(db.Model):
    __tablename__ = 'SNP_risk_stats'
    mapped_gene_stats = db.Column('Mapped_Gene', db.String(100), nullable=True, primary_key=True)
    risk_mean = db.Column('Mean_Risk_Allele_Frequency', db.Float, nullable=True)
    risk_std = db.Column('STD_Risk_Allele_Frequency', db.Float, nullable=True)

# class for tajima D of south asian population
class Tajima(db.Model):
    __tablename__ = 'TajimasD_stats' 
    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column('CHROM', db.Integer, unique=False)
    bin_start = db.Column('BIN_START', db.Integer, unique= False)
    tajD = db.Column('TajimaD', db.Float, nullable=True)
    sa_pop = db.Column('Population', db.String(50))

# class for nSL of south asian population
class NSL(db.Model):
    __tablename__ = 'nSL_stats'
    id = db.Column('rowid', db.Integer, primary_key=True)
    chrom_num = db.Column('Chromosome_Number', db.Integer, unique=False)
    win_start = db.Column('Window_Size_Lower', db.Integer)
    nsl_val = db.Column('Standardised_nSL', db.Float, nullable=True)
    nsl_sa_pop = db.Column('population', db.String(50))
    nsl_pos = db.Column('Physical_Position', db.Integer)
    
# class with population information
class pop_info(db.Model):
    __tablename__ = 'population'
    p_type = db.Column('sampling_location', db.String(50), primary_key=True)
    p_desc = db.Column('population_description', db.String(50))

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
    # default value when the page is first loaded
    snps = None  
    if request.method == 'POST':
        # grab user input
        query1 = request.form.get('query1', '')
        query2 = request.form.get('query2', '')
        query3 = request.form.get('query3', '')
        query4 = request.form.get('query4', '')
            
        # Join 2 tables and query using the input
        # return all information matching any parameter entered
        if query1 or query2 or query3 or query4:
            snps = db.session.query(SNP, SNP_sumstat).outerjoin(SNP_sumstat, SNP.mapped_gene == SNP_sumstat.mapped_gene_stats).filter(
                (SNP.rs_value == query1) | 
                (SNP.gene_pos == query2) | 
                (SNP.snp_population == query4) |
                (SNP.mapped_gene == query3)).all()
                
    return render_template("query.html", snps=snps)


# ontology page
@app.route('/query/ontology/<mapped_gene>/')
def ontology(mapped_gene):
    ontology_results = Ontology.query.filter_by(gene_name=mapped_gene).all()
    # Render the ontology results page
    return render_template('ontology.html', mapped_gene=mapped_gene, ontology_results=ontology_results)

# visualisation page for tajima and nSL
@app.route('/query/visualisation/<rs_value>/', methods=['GET', 'POST'])
def visualisation(rs_value):

    # query into SNP table to grab basic snp information
    snp = db.session.query(SNP).filter_by(rs_value=rs_value).first()
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
        # pop_info_disp = population_info.get(query5, "Please select a population")
        pop_info_disp = pop_info.query.filter_by(p_type=query5).first()
        if pop_info_disp:
             pop_info_disp = pop_info_disp.p_desc 
        else:
            pop_info_disp = "Please select a population"
         
        # if Tajimas or nSL stat is selected...
        if query6 in ["Tajimas_D", "nSL"]:

            # following block of code is if Tajima is selected
            if query6 == "Tajimas_D":
                stat_title = "Tajimas_D"

                # Calculate selection window based off user input and postion
                window = (position // 10000) * 10000
                lower = window - (int(query7) * 1000)
                upper = window + (int(query7) * 1000)
            
                # Query the DB

                # if all populations are selected
                if query5 == "All":
                    filt = db.session.query(Tajima).filter(
                        (Tajima.chrom == chromosome) &
                        (Tajima.bin_start.between(lower, upper))).all()
                # if a specfic population is chosen
                else:
                    filt = db.session.query(Tajima).filter(
                        (Tajima.sa_pop == query5) &
                        (Tajima.chrom == chromosome) &
                        (Tajima.bin_start.between(lower, upper))
                    ).all()

                # once the databse is queried, create the dataframe
                if filt:
                    df = pd.DataFrame([row.__dict__ for row in filt])
                    df.drop(columns=['_sa_instance_state'], inplace=True)

                    # Compute region stats
                    region_mean = round(df['tajD'].mean(), 4)
                    region_std = round(df['tajD'].std(), 4)
                    stat_value = next((row.tajD for row in filt if row.bin_start == window))

                    df["Tajima's D Mean (Selected Region)"] = region_mean
                    df["Tajima's D Std. (Selected Region)"] = region_std  
                    df["Tajima's D"] = stat_value

                    # to show the tajima d value to user
                    final_stat = df.iloc[0,7]

                    # interactive scatter plot

                    # color by population if user selected "All"
                    color_col = 'sa_pop' if query5 == "All" else None
                    # dynamic title
                    if query5 == "All":
                        plot_title = f"Tajima's D for Chromosome Window {window} ± {query7} kb (All Populations)"
                    else:
                        plot_title = f"Tajima's D for Chromosome Window {window} ± {query7} kb ({query5} Population)"

                    # plot based on df
                    fig = px.scatter(
                        df,
                        x = 'bin_start',
                        y = 'tajD',
                        color = color_col,  # color by population only if "All"
                        title = plot_title,
                        hover_data = ["bin_start", "tajD"]
                    )
                    # update the hover data to use more readable
                    fig.update_traces(
                        hovertemplate=(
                        "<b>Chromosome Position:</b> %{x}<br>"  # Bin start (Chromosome Position)
                        "<b>Tajima's D Value:</b> %{y}<br>"  # Tajima's D value
                    )
                    )
                    # update legend title if all populations selected
                    if query5 == "All":
                        fig.update_layout(
                            legend_title="Populations:"  # Change legend title
                        )
                    fig.update_layout(
                        xaxis_title = f"Chromosome {chromosome} Region (bp)",
                        yaxis_title = "Tajima's D"
                    )

                    # add lines (threshold, region mean)
                    fig.add_hline(y=-2, line_color="red", annotation_text="Threshold (-2)", annotation_position="bottom left")
                    fig.add_hline(y=region_mean, line_color="green", annotation_text="Region Mean", annotation_position="bottom left")

                    # convert to HTML for template rendering
                    plot_html = fig.to_html(full_html=False)

                    # handle the Download request
                    if 'download' in request.form:
                        tsv_buffer = BytesIO()
                        df.to_csv(tsv_buffer, index=False, sep="\t")
                        tsv_buffer.seek(0)
                        return send_file(
                            tsv_buffer,
                            mimetype="text/tab-separated-values",
                            as_attachment=True,
                            download_name=f"Tajimas_D_{rs_value}.tsv"
                        )

                    # return the plot HTML in the template
                    return render_template(
                        'visualisation.html',
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


            # following block of code is if Tajima is selected
            elif query6 == "nSL":
                stat_title = "nSL"

                # Calculate window boundaries
                lower = position - (int(query7) * 1000) 
                upper = position + (int(query7) * 1000)
                window = position

                # Query the DB
                if query5 == "All":
                    filt = db.session.query(NSL).filter(
                        (NSL.chrom_num == chromosome) &
                        (NSL.nsl_pos.between(lower, upper))
                    ).all()
                else:
                    filt = db.session.query(NSL).filter(
                        (NSL.nsl_sa_pop  == query5) &
                        (NSL.chrom_num == chromosome) &
                        (NSL.nsl_pos.between(lower, upper))
                    ).all()

                if filt:
                    df = pd.DataFrame([row.__dict__ for row in filt])
                    df.drop(columns=['_sa_instance_state'], inplace=True)

                    # compute region stats
                    df["search_position"] = position
                    region_mean = round(df['nsl_val'].mean(), 4)
                    region_std = round(df['nsl_val'].std(), 4)
                    stat_value = next(
                        (getattr(row, 'nsl_val') 
                        for row in filt 
                        if getattr(row, 'win_start') == window),
                        None
                    )
                    # then add to dataframe in new columns
                    df["nSL Mean (Selected Region)"] = region_mean
                    df["nSL Std. (Selected Region)"] = region_std  
                    df["nSL"] = stat_value

                    # value to present to user as the nSL value on SNP of interest
                    final_stat = final_stat = df.loc[df['nsl_pos'] == position, 'nsl_val'].iloc[0] if not df.loc[df['nsl_pos'] == position, 'nsl_val'].empty else None

                    # interactive scatter plot

                    # color by population if user selected "All"
                    color_col = 'nsl_sa_pop' if query5 == "All" else None
                    # dynamic title
                    if query5 == "All":
                        plot_title = f"nSL for Chromosome Window {window} ± {query7} kb (All Populations)"
                    else:
                        plot_title = f"nSL for Chromosome Window {window} ± {query7} kb ({query5} Population)"


                    # plot based on dataframe from query
                    fig = px.scatter(
                        df,
                        x='nsl_pos',
                        y='nsl_val',
                        color=color_col,  # color by population only if "All"
                        title=plot_title,
                        hover_data=["win_start", "nsl_val"]
                    )
                    fig.update_layout(
                        xaxis_title=f"Chromosome {chromosome} Region (bp)",
                        yaxis_title='nSL',
                        xaxis=dict(
                        range=[lower, upper],
                        showgrid=True,
                        zeroline=True,
                        )
                    )
                    
                    # update legend title if all populations select
                    if query5 == "All":
                        fig.update_layout(
                            legend_title="Populations:"  # Change legend title
                        )
                    # update the hover data to use more readable
                    fig.update_traces(
                        hovertemplate=(
                        "<b>Chromosome Position:</b> %{x}<br>"  # Chromosome Position
                        "<b>nSL Value:</b> %{y}<br>"  # nSL value
                        )
                    )
                    # add horizontal lines (threshold, region mean)
                    fig.add_hline(y=-2, line_color="red", annotation_text="Threshold (+2)", annotation_position="bottom left")
                    fig.add_hline(y=2, line_color="red", annotation_text="Threshold (-2)", annotation_position="bottom left")
                    fig.add_hline(y=region_mean, line_color="green", annotation_text="Region Mean", annotation_position="bottom left")
                    # vertical line for snp of interest
                    fig.add_vline(x=position, line=dict(color="black", width=2, dash="dash"))
                    fig.add_annotation(
                        x=position, y=1,  
                        xref="x", yref="paper",  # y coordinate in paper space (0 to 1)
                        text="Position of SNP chosen",
                        showarrow=True, arrowhead=2,
                        ax=0, ay=-40,      
                        font=dict(size=12, color="black"),
                        align="center", yanchor="top"  
                    )


                    # convert to HTML for template rendering
                    plot_html = fig.to_html(full_html=False)

                    # handle the Download request
                    if 'download' in request.form:
                        tsv_buffer = BytesIO()
                        df.to_csv(tsv_buffer, index=False, sep="\t")
                        tsv_buffer.seek(0)
                        return send_file(
                            tsv_buffer,
                            mimetype="text/tab-separated-values",
                            as_attachment=True,
                            download_name=f"nSL_{rs_value}.tsv"
                        )

                    # return the plot HTML in the template
                    return render_template(
                        'visualisation.html',  
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

    # GET request
    return render_template('visualisation.html', rs_value=rs_value, snp=snp)


# for downloading the data
@app.route('/query/visualisation/<rs_value>/download', methods=['GET'])
def download_stats(rs_value):

    # grab the relevant info from url
    stat_type = request.args.get('stat', 'Tajimas_D')  # "Tajima" or "nSL"
    pop = request.args.get('pop', 'All') # population info
    window_size = int(request.args.get('window', 10)) # default is 10

    snp = db.session.query(SNP).filter_by(rs_value=rs_value).first()
    if not snp:
        return "SNP not found", 404 # error handling
    chromosome = snp.chr_id
    position = snp.gene_pos

    lower = position - (window_size * 1000)
    upper = position + (window_size * 1000)

    # code below for taj & nSL is duplicate of previously defined functions
    if stat_type == "Tajimas_D":
        query = db.session.query(Tajima).filter(
            Tajima.chrom == chromosome,
            Tajima.bin_start.between(lower, upper)
        )
        if pop != "All":
            query = query.filter(Tajima.sa_pop == pop)

        filt = query.all()
        if not filt:
            return "No Tajima's D data found.", 404

        df = pd.DataFrame([row.__dict__ for row in filt])
        df.drop(columns=['_sa_instance_state'], inplace=True)

        df["Tajima's D Mean (Region)"] = df['tajD'].mean()
        df["Tajima's D Std (Region)"] = df['tajD'].std()

        filename = f"Tajimas_D_{rs_value}.tsv"

    elif stat_type == "nSL":
        query = db.session.query(NSL).filter(
            NSL.chrom_num == chromosome,
            NSL.nsl_pos.between(lower, upper)
        )
        if pop != "All":
            query = query.filter(NSL.nsl_sa_pop == pop)

        filt = query.all()
        if not filt:
            return "No nSL data found.", 404

        df = pd.DataFrame([row.__dict__ for row in filt])
        df.drop(columns=['_sa_instance_state'], inplace=True)

        df["nSL Mean (Region)"] = df['nsl_val'].mean()
        df["nSL Std (Region)"] = df['nsl_val'].std()

        filename = f"nSL_{rs_value}.tsv"

    else:
        return "Unknown statistic type.", 400  # error handling

    # convert df to tsv file
    tsv_buffer = BytesIO()
    df.to_csv(tsv_buffer, index=False, sep="\t")
    tsv_buffer.seek(0)

    # downloads the file for user
    return send_file(
        tsv_buffer,
        mimetype="text/tab-separated-values",
        as_attachment=True,
        download_name=filename

    )

def fetch_snp_sequence(rsid):
    """Fetch SNP location, sequence, and modify the SNP position using mapping keys."""
    base_url = f"https://rest.ensembl.org/variation/homo_sapiens/{rsid}?content-type=application/json"
    response = requests.get(base_url)
    if response.status_code == 200:
        data = response.json()
        mapping = data["mappings"][0]
        allele_string = mapping.get("allele_string", "N/A")
        print("HELLO!!!!")
        allele_string = mapping.get("allele_string", "N/A")
        print(allele_string)
        # Check if mappings exist and have data
        if "mappings" not in data or len(data["mappings"]) == 0:
            print(f"No SNP mapping data found for {rsid}")
            return None

        # Use the first mapping for the SNP info
        mapping = data["mappings"][0]
        if "seq_region_name" not in mapping or "start" not in mapping:
            print(f"No SNP data found for {rsid} in mapping")
            return None

        chrom = mapping["seq_region_name"]
        snp_position = mapping["start"]

        # Set flank size to 500 on either side of the SNP
        flank_size = 500
        region_start = max(1, snp_position - flank_size)
        region_end = snp_position + flank_size

        sequence_url = (
            f"https://rest.ensembl.org/sequence/region/human/"
            f"{chrom}:{region_start}..{region_end}?content-type=text/plain"
        )
        seq_response = requests.get(sequence_url)
        if seq_response.status_code == 200:
            sequence = seq_response.text.strip()
            result = {
                "rsid": rsid,
                "chromosome": chrom,
                "start": region_start,
                "snp_position": snp_position,
                "allele_string" : allele_string,
                "sequence": sequence
                 
            }
            return result
        else:
            print(f"Failed to fetch sequence for {rsid}")
    else:
        print(f"Failed to fetch SNP info for {rsid}")
    return None

@app.route('/query/visualisation/<rs_value>/sequencevisualisation', methods=['GET'])
def sequence_visualisation(rs_value):
    snp_data = fetch_snp_sequence(rs_value)
    if snp_data:
        print("hello!!")
        print(json.dumps(snp_data, indent=6))  # Pretty print to check structure
        return render_template(
            "sequence_visualisation.html",
            snp_data=snp_data,
            region_start=snp_data["start"],
            snp_position=snp_data["snp_position"]
            
        )
    
    
    return "Failed to fetch SNP data", 500

if __name__ == '__main__':
    app.run(debug=True)
