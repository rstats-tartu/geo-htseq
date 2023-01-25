[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7529832.svg)](https://doi.org/10.5281/zenodo.7529832)


# README #


### What is this repository for? ###

Assessing the quality of experiments, both directly and through replication attempts, is becoming a grave concern in biomedicine. This is especially true in omics experiments, where thousands of independent measurements are done in parallel, with the understanding that only a small minority will probe scientifically meaningful effects – a state of affairs conductive to mistaking false positives for scientific discoveries.  

Here, we assessed the quality of the high-throughput sequencing experiments submitted to Entrez GEO database until 2020-12-31.

### Descripton of the workflow ###

NCBI GEO datasets queries were performed using Bio.Entrez python package and by sending requests to NCBI Entrez public API. The query was 'expression profiling by high throughput sequencing[DataSet Type] AND ("2000-01-01"[PDAT] : "2019-12-31"[PDAT])'. FTP links from GEO datasets document summaries were used to download supplementary files lists. Supplementary files were filtered for downloading, based on file extensions, to keep file names with "tab", "xlsx", "diff", "tsv", "xls", "csv", "txt", "rtf", and "tar" file extensions. We dropped files whos names contained regular expression "filelist.txt|raw.tar$|readme|csfasta|(big)?wig|bed(graph)?|(broad_)?lincs", because we were not expecting to find P values from these files. Downloaded files were imported using python pandas package, and searched for unadjusted P value sets. Unadjusted P value sets and summarised expression level of associated genomic features were identified using column names. P value columns from imported tables were identified by using regular expression "p[^a-zA-Z]{0,4}val", adjusted P value sets were identified using regular expression "adj|fdr|corr|thresh" and omitted from furter analysis. Expression levels of genomic features were identified by using following regular expressions: "basemean", "value", "fpkm", "logcpm", "rpkm", "aveexpr". When multiple expression level columns were present in a table, then average expression level was calculated for each feature. 
Identified raw P value sets were classified based on their histogram shape. 
Histogram shape was determined based on the presence and location of peaks.
P value histogram peaks (bins) were detected using a quality control threshold described in [1], a Bonferroni-corrected $\alpha$-level quantile of the cumulative function of the binomial distribution with size m and probability p. Histograms, where none of the bins were over QC-threshold, were classified as "uniform". Histograms, where bins over QC-threshold started either from left or right boundary and did not exceeded 1/3 of the 0 to 1 range, were classified as "anti-conservative" or "conservative", respectively. Histograms with peaks or bumps in the middle or with non-continuous left- or right-side peaks were classified as "other". Histograms with peaks on both left- and right-side were classified as "bimodal".
Raw P value sets with anti-conservative shape were used to calculate the pi0 statistic using smoother method [2] as implemented in gdsctools Python library. When expression level data were identifiable in imported tables, raw P values were further filtered to remove uninformative features using following thresholds: basemean=10, logcpm=1, rpkm=1, fpkm=1, aveexpr=3.32. Differential expression analysis platform were inferred from column name pattern for cuffdiff, DESeq/DESeq2, EdgeR, and limma, all other unidentified sets were binned as "unknown". Publication data were downloaded from NCBI PubMed database. Citation data were downloaded from Elevier Scopus database. Sequence read library statistics were downloaded from NCBI SRA database. The code is available as a snakemake workflow on [tpall/geo-htseq](https://github.com/tpall/geo-htseq) Github repo. All versions of the dataset is deposited in Zenodo with doi:10.5281/zenodo.3859722.

### How do I get set up? ###

* To get started you need to download and install miniconda3 and create conda environment with snakemake
* Go to <https://docs.conda.io/en/latest/miniconda.html> for download and installation instructions of miniconda3
* Create conda environment with snakemake, essentially as described in <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>

```bash
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

* (Fork and) clone this repository 

```bash
git clone https://github.com/rstats-tartu/geo-htseq.git
```

* cd to working directory and activate conda environment

```bash
cd geo-htseq
conda activate snakemake
```

* Dry run workflow

```bash
snakemake -n
```

* To run the workflow you need to set up
    - [NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) as NCBI_APIKEY environment variable
    - [Elsevier API key](https://dev.elsevier.com) as ELSEVIER_GEOSEQ environment variable


* Run workflow. 

```bash
snakemake --use-conda -j
```


### Who do I talk to? ###

* Taavi Päll taavi.pall@ut.ee

### References ###
[1] Breheny, P., Stromberg, A., & Lambert, J. (2018). p-Value Histograms: Inference and Diagnostics. High-throughput, 7(3), 23. https://doi.org/10.3390/ht7030023    
[2] J. D. Storey and R. Tibshirani. Statistical significance for genome-wide experiments. Proceedings of the National Academy of Sciences, 100:9440–9445, 2003. https://www.pnas.org/content/100/16/9440
