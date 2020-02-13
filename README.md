# README #


### What is this repository for? ###

Assessing the quality of experiments, both directly and through replication attempts, is becoming a grave concern in biomedicine. This is especially true in omics experiments, where thousands of independent measurements are done in parallel, with the understanding that only a small minority will probe scientifically meaningful effects – a state of affairs conductive to mistaking false positives for scientific discoveries.  

Here, we assessed the quality of the high-throughput sequencing experiments submitted to Entrez GEO database until 2018-12-31.

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