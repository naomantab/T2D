# Project Overview

This software was developed by "Jammy Jelly" student group of Queen Mary University of London for the purpose of a software development projects as part of the MSc Bioinformatics Program.

This project involves the development of a **web-based software tool** to handle **molecular biology data**. The primary goal is to create an application that retrieves and integrates **single nucleotide polymorphisms (SNPs)** associated with **Type 2 Diabetes (T2D)** and visualizes relevant population genomic and functional data.

## 🧪 Website features
- Retrieve **SNP information** based on:
  - SNP ID (rs value)
  - Genomic coordinates (chromosome, start, and end)
  - Mapped gene name
  - Population
    
- Display the following details for each SNP:
  - SNP name
  - genomic position
  - p-value from genome-wide association study
  - Mapped gene(s)
  - Summary statistics of **positive selection** for at least two South Asian populations
- **Gene ontology analysis:**
  - Select a mapped gene to retrieve functional/ontology terms
- **Population-level analysis:**
  - Select populations of interest
  - Display **descriptive statistics**
  - **Visualize summary statistics** for positive selection across genomic regions
- **Exportable results:**
  - Download a text file (TSV) containing summary statistics (average & standard deviation)

## 🖥️ Technologies Used
- **Frontend:** HTML, CSS, JavaScript
- **Backend:** Python (Flask)
- **Database:** SQL (SQLite)
  https://drive.google.com/drive/folders/1J-4i_U8hxD_Z2teKKpeK9lQdk0ca0Ugc?usp=sharing
- **Data Sources:**
  - [GWAS Catalog](https://www.ebi.ac.uk/gwas)
  - [T2D Knowledge Portal](https://t2d.hugeamp.org)
  - [International Genome Sample Resource](https://www.internationalgenome.org/)
  - [Gene Ontology Database](https://geneontology.org/)
- **Visualization:** Plotly


## ⚙️ Installation & Usage
    Download Database:
    ```sh
    https://bit.ly/t2d-db
    ```
1. Clone the repository:
   ```sh
   git clone https://github.com/naomantab/T2D.git
   cd T2D
   ```
2. Set up a virtual environment:
   ```sh
   python -m venv env
   source env/bin/activate  # On Windows use: env\Scripts\activate
   ```
3. Run the application:
   ```sh
   flask run  # For Flask
   ```
## ⭐ Project Team & Acknowledgments

This project is part of the **MSc Bioinformatics Software Development Group Project 2025**.  
Special thanks to **Professor Conrad Bessant** and **Dr Matteo Fumagalli** for their guidance.

**Contributors:**
- Naoman Tabassam
- Deepak Raj
- Katelin Cunningham
- Olayemi Bakare
- Joseph Eytle


## ⚠️ License
N/A for now

   
