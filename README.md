# E2A Paper
Code for "Interpretable deep learning reveals the role of an E-box motif in suppressing somatic hypermutation of AGCT motifs within human immunoglobulin variable regions" by Tambe, MacCarthy and Pavri published in Frontiers in Immunology.

## Usage
Clone the repository
```bash
git clone https://github.com/abhikt/e2a_paper
```
Run the pipeline script to generate data file for plotting
```bash
python run_pipeline.py
```
Plot the figures using
```bash
python plot_figures.py
```
Make sure to unzip the large data files
```bash
gunzip data/gm12878_e2achip_rpkm_500.csv.gz
gunzip data/ramos_e2achip_rpkm_500.csv.gz
```
