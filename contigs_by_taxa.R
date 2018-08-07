## Select contigs based on taxonomic assignments
## Chunyu Zhao 2018/05/22

library(dplyr)
library(readr)
if (! "taxonomizr" %in% rownames(installed.packages()))
  install.packages("taxonomizr")
library(taxonomizr)

## input
taxaNodes <- read.nodes(file.path(snakemake@input[['taxa_nodes']]))
taxaNames <- read.names(file.path(snakemake@input[['taxa_names']]))
accessionTaxasql <- file.path(snakemake@input[['taxa_sql']])

## parameters
db <- snakemake@params[['db']]
taxaName <- snakemake@params[['taxaName']]
min_contig_len <- as.numeric(snakemake@params[['min_contig_len']])

## output
summary_files <- snakemake@input[['summary_files']]
outdir <- snakemake@params[['outdir']]
report_fp <- snakemake@output[['report']]

## parse blastn  
blastn <- do.call(rbind, 
                  lapply(1:length(summary_files), 
                         function(x) read_delim(summary_files[x], delim="\t"))) %>% 
  filter(length >= min_contig_len) %>%
  filter_(paste("!is.na(", db, ")", sep="")) %>%
  mutate(taxaID = accessionToTaxa(get(db), accessionTaxasql)) 

blastn <- cbind(blastn, getTaxonomy(blastn$taxaID,taxaNodes,taxaNames))

## filter by taxaName, this is slower, yet RefSeq doesn't necessarily update the new assembly 
taxaName <- gsub("_", " ", taxaName)
blastn_taxa <- filter(blastn, grepl(taxaName, species, perl=T, ignore.case = T))

## write the contig names of the filtered contigs into file per sample

dir.create(outdir, showWarnings = FALSE)
customWrite <- function(DF, outdir=".") {
  fname <- file.path(outdir, paste(unique(DF$sample),".txt", sep=""))
  print(fname)
  write.table(DF$contig, fname, row.names = F, col.names = F, quote = F)
  return(DF)
}

done <- blastn_taxa %>%
  group_by(sample) %>%
  do(customWrite(., outdir=outdir))

## Caution: not all samples will have the taxa of interest
write.table(done, report_fp, sep="\t", row.names=F, quote=F)
