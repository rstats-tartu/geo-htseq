# Filter downloaded supplementary file names ------------------------------
# In this section, we filter supplementary file names for patterns: 
# we are looking only for tabular data. 
out_string1 <- c("filelist","annotation","readme","error","raw.tar","csfasta",
                 "bam","sam","bed","[:punct:]hic","hdf5","bismark","map",
                 "barcode","peaks")
out_string2 <- c("tar","gtf","(big)?bed(\\.txt|12|graph|pk)?","bw","wig",
                 "hic","gct(x)?","tdf","gff(3)?","pdf","png","zip","sif",
                 "narrowpeak","fa(sta|stq)?", "r$", "rda(ta)?$")