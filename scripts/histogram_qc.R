library(ggplot2)

histogram_qc <- function(x, breaks, fdr = 0.05) {
  
  x <- na.omit(x)
  b <- 1 / breaks
  h <- hist(x, breaks = seq(0, 1, b), plot = FALSE)
  counts <- h$counts
  mids <- h$mids
  qc <- qbinom(1 - b * fdr, length(x), b)
  counts_over_qc <- counts > qc
  mids_over_qc <- mids[counts_over_qc]
  
  if (all(!counts_over_qc)) {
    
    Class <- "uniform"
    QC_pass <- TRUE
    
  } else {
    
    end_of_peak <- rle(counts_over_qc)$lengths[1] + 2
    QC_pass <- !any(counts_over_qc[end_of_peak:length(counts_over_qc)])
    
    if (QC_pass) {
      
      Class <- "anti-conservative"
      
    } else {
      
      end_of_peak <- rle(rev(counts_over_qc))$lengths[1] + 2
      conservative <- !any(rev(counts_over_qc)[end_of_peak:length(counts_over_qc)])
      
      if (conservative) {
        
        Class <- "conservative"
        
      } else {
        
        Class <- "other"
      }
    }
  }
  
  list(QC_pass, Class, mids_over_qc)
}

qc_plot <- function(x, breaks) {
  x <- na.omit(x)
  b <- 1 / breaks
  qc <- qbinom(1 - b * 0.05, length(x), b)
  ggplot(data = NULL) +
    geom_histogram(aes(x), binwidth = 1 / breaks, center = 1 / (2 * breaks)) +
    geom_hline(yintercept = qc, linetype = "dashed") +
    theme(axis.title = element_blank())
}

