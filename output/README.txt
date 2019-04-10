
R/results_pvalues.R script generates following csv files: 
- output/taxons.csv -- taxonomic metadata for GEO series 
   - Accession -- GEO accession number
   - organism -- name of the organism
   - taxid -- taxon id of the organism
   - multiple taxons -- serie includes also other taxa, only first is shown

- output/number_of_samples_in_geo.csv -- tally GEO series with 1..n samples
   - samples -- number of samples in GEO serie
   - N -- count of GEO series

- output/pvalues_pool_pub.csv -- set of vars from pvalues_pool data frame for use in results_publication.R
   - Accession -- GEO accession number
   - suppdata_id -- GEO supplementary table id
   - Type -- class assigned to pvalue histogram
   - pi0 -- proportion of true nulls calculated by limma::propTrueNull()  
   - SRP -- SRP calculated by srp function
   - pi01 -- proportion of true nulls calculated by qvalue package   
   - fp -- number of false positive effects
   - rs -- don't remember..
   - ud -- number of undetected effects

output/pvalue_sets_classes.csv -- table for manual reclassification.

output/srp_stats.csv -- srp stats for anti-conservative pvalue histograms for raw and basemean filtered sets

R/pvalue_classification.R script generates following csv file, imported by R/results_pvalues.R: 
- output/pvalue_histogram_nnet_classification.csv -- nnet classification results
   - suppdata_id -- GEO supplementary table id
   - Filter -- raw or basemean filtered pvalue sets flag 
   - human  -- histogram class assigned by Ülo
   - nnet -- histogram class assigned by nnet model 
   

