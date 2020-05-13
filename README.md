# README #


### What is this repository for? ###

Assessing the quality of experiments, both directly and through replication attempts, is becoming a grave concern in biomedicine. This is especially true in omics experiments, where thousands of independent measurements are done in parallel, with the understanding that only a small minority will probe scientifically meaningful effects – a state of affairs conductive to mistaking false positives for scientific discoveries.  

Here, we assessed the quality of the high-throughput sequencing experiments submitted to Entrez GEO database until 2019-12-31.

### Descripton of the workflow ###

NCBI GEO datasets queries were performed using Bio.Entrez python package and by sending requests to NCBI Entrez public API. FTP links from GEO datasets document summaries were used to download supplementary files lists. Supplementary files were filtered for downloading, based on file extensions, to keep file names with "tab", "xlsx", "diff", "tsv", "xls", "csv", "txt", "rtf", and "tar" file extensions. Downloaded files were imported using python pandas package, and searched for unadjusted P value sets. Unadjusted P value sets and summarised expression level of associated genomic features were identified using column names. Identified raw P value sets were first classified based on their histogram shape. Raw P value sets with anti-conservative shape were used to calculate pi0 statistic. When expression level data were present and identifiable in imported tables, raw P values were further filtered to remove uninformative features using following thresholds basemean=10, logcpm=1, rpkm=1, fpkm=1, aveexpr=3.32. Differential expression analysis platform were inferred using column name pattern for cuffdiff, DESeq/DESeq2, EdgeR, and limma, all other unidentified sets were binned as unknown/unidentified. Publication data were downloaded from NCBI PubMed database. Citation data were downloaded from Elevier Scopus database. Sequence read library statistics were downloaded from NCBI SRA database. The code is available as a snakemake workflow on [tpall/geo-htseq](https://github.com/tpall/geo-htseq) Github repo. Dataset is deposited in Zenodo with doi:10.5281/zenodo.3778160.

### How do I get set up? ###

* To get started you need to download and install miniconda3 and create conda environment with snakemake
* Go to <https://docs.conda.io/en/latest/miniconda.html> for download and installation instructions of miniconda3
* Create conda environment with snakemake

```bash
conda create -n snakemake-env -c bioconda snakemake
```

* (Fork and) clone this repository 

```bash
git clone https://github.com/tpall/geo-htseq.git
```

* Dry run workflow

```bash
cd geo-rnaseq
conda activate snakemake-env
snakemake -n
```

* To run the workflow you need to set up
    - [NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) as NCBI_APIKEY environment variable
    - [Elsevier API key](https://dev.elsevier.com) as ELSEVIER_GEOSEQ environment variable


* We have been running this workflow in SLURM cluster using following command. *cluster.json* file contains cluster configuration info (e.g. time, mem) for jobs. 

```bash
snakemake --use-conda --cluster-config cluster.json --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --output {cluster.output}" -j
```


### Who do I talk to? ###

* Taavi Päll taavi.pall@ut.ee