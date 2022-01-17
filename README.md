# cryptic_donor_prediction
code to generate cryptic donor database figures in Dawes et al. 2022 (Empirical prediction of variant-activated cryptic splice donors using population-based RNA-Seq data)

# Overview
- `data/` contains raw and reference files required to run analyses

- The cryptic-donor database is available pre-computed `data/processed/css_database.tsv.gz`

- `src/processing` contains R notebooks used to process data

- `notebook/` contains R notebooks (`figures_code/`) to generate figures in the paper, as well as the outputted figures (`figures/`). Intermediate processed files are required to run these notebooks - running data processing notebooks in `src/processing` is therefore a required precursor.


### Data processing notebooks

1. cryptic variant processing:
`src/processing/vd_cryptic_processing.Rmd` processes the source variants located in `data/raw/combined_cryptic_variants.tsv`, adding surrounding sequence and finding REF and ALT authentic, cryptic and decoy-donor sequences, and appends predicted strength scores from 4 algorithms (MES, NNS, SAI and DF)
Output : `data/processed/css_database.tsv.gz` 
`data/` contains a legend for the css database file.

2. cryptic variant splice competent events processing:
`src/processing/cssDb_missplicing.Rmd` takes the cryptic database (`data/processed/css_database.tsv.gz`) and appends information about events found in 40K-RNA
Output: `data/processed/css_database_msd.tsv.gz`

3. cryptic variant algorithm predictions processing:
`src/processing/cssDb_predictions.Rmd` takes the cryptic database as processed in step 2 (`data/processed/css_database_msd.tsv.gz`) and categorises the predictions given by each of the four algorithms for use in Figure 3
Output: `data/processed/css_database_predictions.tsv.gz`

4. Human genome decoys processing:
`src/processing/hg19_processing.Rmd` processes the annotated human donors in `data/ref/ensembl75_canonical_[introns/exons].tsv` searching within up to 500nt of surrounding sequence for all possible decoy donors and annotating with DF scores.
Output: `data/processed/human_introns_decoys.tsv.gz` 
`data/` contains a legend for the human introns decoys file.

5. Human genome decoys simulations processing: 
`src/processing/hg19_seq_shuffle_decoys.Rmd` performs the process described in Figure S4 of Dawes et al. 2022, shuffling exonic and intronic sequences and tallying decoys across 15 simulations.
Output: `data/processed/human_introns_decoys_simulations.tsv.gz`

### Figures

Once the processing notebooks have been run, the code to create figures (`notebook/figures_code/`) can be run and figures can be produced. 


