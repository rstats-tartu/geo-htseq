
n <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_size <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX"))

if (is.na(n) & is.na(array_size)) {
  array_size <- n <- 1
}

# Load libs, settings and functions
source("R/_common.R")
source("R/munge_geo.R")
source("R/checkFullRank.R")
source("R/text_funs.R")
pacman::p_load_gh("seandavi/GEOquery")

munge_suppfiles <- function(docsums, gsem, n, array_size, out_path) {
  # GEO query results and document summaries -----
  ds <- read_rds(docsums)
  
  # Load series matrix data frames -----
  gsem <- read_rds(gsem)
  
  # Extract series matrixes
  gsem <- gsem %>% 
    mutate(series_matrix = map(gse, "result"))
  
  # Identify and filter out Accessions with missing gsematrices
  gsem_missing_or_faulty <- gsem %>% 
    filter(map_lgl(series_matrix, ~class(.x) != "ExpressionSet"))
  
  gsem <- gsem %>% 
    filter(map_lgl(series_matrix, ~class(.x) == "ExpressionSet"))
  
  ## Read in local supplemental tables -----
  
  # File name extensions in downloaded supplementary files
  supptabs <- geofile_df(local_suppfile_folder, "suppfiles")
  
  # Unzip gz xls(x)? files, keep originals
  xls_gz <- filter(supptabs, str_detect(suppfiles, "xls(x)?.gz$"))$suppfiles
  xls_gz <- file.path(local_suppfile_folder, xls_gz)
  
  # List files again to exclude compressed xls files
  supptabs <- geofile_df(local_suppfile_folder, "suppfiles") %>% 
    filter(!str_detect(suppfiles, "xls(x)?.gz$"))
  
  # Duplicated files
  supptabs_duplicated <- supptabs %>% 
    mutate(files = str_replace(suppfiles, "\\.(gz|bz2)$","")) %>% 
    filter(duplicated(files))
  
  # Merge gsem ExpressionSets to supptabs
  supptabs <- select(gsem, Accession, series_matrix) %>% 
    nest(series_matrix, .key = "matrixfiles") %>% 
    left_join(supptabs, .)
  
  # Files that crash R session
  bad <- c("GSE93374_Merged_all_020816_DGE.txt.gz", 
           "GSE88931_RNA-seq_MergedReadCounts.tsv.gz",
           "GSE55385_transcripts_GSE.tsv.gz",
           "GSE74549_ChIP_1kbWindows_correctedReadCount.txt.gz",
           "GSE77213_Nguyen_GEO_TN03_tallies.xls",
           "GSE77213_Nguyen_GEO_TN05_tallies_total.xls",
           "GSE60012_100bpTiles_RRBS_Mouse.txt.gz",
           "GSE60415_heatmap-upload.xls",
           "GSE67516_RNA_seq_rep1_diffExp_analysis.xls",
           "GSE53298_processed_data.xlsx",
           "GSE53350_SAGE_KomatsuFurukawa_Processed.xls",
           "GSE60483_AQ_exp.tab.gz",
           "GSE53260_Supplementary_S2.xls",
           "GSE53260_Supplementary_S1.xls",
           "GSE89113_T35_vs_T45.T45.UP.0.05.A3SS.MATS.JunctionCountOnly.txt.gz",
           "GSE89113_T35_vs_T45.T45.UP.0.05.MXE.MATS.JunctionCountOnly.txt.gz",
           "GSE89113_T40_vs_T45.A3SS.MATS.JunctionCountOnly.txt.gz",            
           "GSE89113_T40_vs_T45.T40.UP.0.05.A3SS.MATS.JunctionCountOnly.txt.gz",
           "GSE89113_T40_vs_T45.T45.UP.0.05.A3SS.MATS.JunctionCountOnly.txt.gz",
           "GSE89113_T35_vs_T45.A3SS.MATS.JunctionCountOnly.txt.gz",            
           "GSE89113_T35_vs_T45.MXE.MATS.JunctionCountOnly.txt.gz",             
           "GSE89113_T25_vs_T40.A3SS.MATS.JunctionCountOnly.txt.gz",            
           "GSE99484_Sulfolobus_acidocaldarius_DSM_639_CP000077_.gbk.txt.gz",   
           "GSE89113_T25_vs_T40.T25.UP.0.05.A3SS.MATS.JunctionCountOnly.txt.gz",
           "GSE89113_T25_vs_T40.T40.UP.0.05.A5SS.MATS.JunctionCountOnly.txt.gz",
           "GSE89113_T25_vs_T45.A3SS.MATS.JunctionCountOnly.txt.gz",            
           "GSE89113_T25_vs_T45.MXE.MATS.JunctionCountOnly.txt.gz",             
           "GSE89113_T25_vs_T45.SE.MATS.JunctionCountOnly.txt.gz",              
           "GSE89113_T25_vs_T45.T25.UP.0.05.SE.MATS.JunctionCountOnly.txt.gz",  
           "GSE89113_T25_vs_T45.T45.UP.0.05.MXE.MATS.JunctionCountOnly.txt.gz", 
           "GSE89113_T35_vs_T40.A5SS.MATS.JunctionCountOnly.txt.gz",            
           "GSE89113_T35_vs_T40.T40.UP.0.05.A5SS.MATS.JunctionCountOnly.txt.gz")
  
  # Remove 'bad' files
  supptabs <- supptabs %>% filter(!(suppfiles %in% bad))
  
  split_supptabs <- split(supptabs, 1:array_size)
  
  split_supptabs[[n]] %>%
    mutate(result = map(suppfiles, ~ try(munge_geo_pvalue(file.path(local_suppfile_folder, .x))))) %>% 
    write_rds(path = glue(out_path))
}

munge_suppfiles(snakemake@input[["docsums"]], 
                snakemake@input[["gsem"]], 
                n = n, 
                array_size = array_size, 
                out_path = snakemake@output[[1]])
