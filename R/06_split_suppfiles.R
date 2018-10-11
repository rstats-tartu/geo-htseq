
# Load libs, local_suppfile_folder, and helper function geofile_df
source("R/_common.R")

munge_suppfiles <- function(docsums, suppfilenames_filtered, gsem) {
  
  # GEO query results and document summaries -----
  ds <- read_rds(docsums)
  
  # supplementary files of interest expected to be locally available
  suppfilenames_filtered <- read_rds(suppfilenames_filtered)
  
  # Load series matrix data frames -----
  gsem <- read_rds(gsem)
  
  # Extract series matrixes
  gsem <- mutate(gsem, series_matrix = map(gse, "result"))
  
  # Filter Accessions with gsematrices present
  gsem <- filter(gsem, map_lgl(series_matrix, ~class(.x) == "ExpressionSet"))
  
  ## Read in local supplemental tables -----
  
  # File names of downloaded supplementary files
  supptabs <- geofile_df(local_suppfile_folder, "suppfiles")
  
  # Keep only files present in filtered suppfilenames
  supptabs <- filter(suppfilenames_filtered, type == "suppl") %>% 
    rename(suppfiles = files) %>% 
    inner_join(supptabs)
  
  # Merge gsem ExpressionSets to supptabs
  supptabs <- dplyr::select(gsem, Accession, series_matrix) %>% 
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
  supptabs <- filter(supptabs, !(suppfiles %in% bad))
  
  mutate(supptabs, splits = rep_len(1:100, nrow(supptabs))) %>% 
    group_by(splits) %>% 
    nest() %>% 
    mutate(data = map2(data, splits, ~write_rds(.x, path = glue::glue("output/tmp/supptabs_{.y}.rds"))))
}

munge_suppfiles(docsums = snakemake@input[["docsums"]], 
                suppfiles_filtered = snakemake@input[["suppfilenames_filtered"]],
                gsem = snakemake@input[["gsem"]])
