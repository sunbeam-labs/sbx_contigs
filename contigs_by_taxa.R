## Select contigs based on taxonomic assignments

library(tidyverse)
library(taxonomizr)

## Inputs:
taxaNodes <- read.nodes(file.path(snakemake@input[['taxa_nodes']]))
taxaNames <- read.names(file.path(snakemake@input[['taxa_names']]))
accessionTaxasql <- file.path(snakemake@input[['taxa_sql']])

db <- snakemake@params[['db']]
taxaName <- snakemake@params[['taxaName']]

summary_files <- snakemake@input[['summary_files']]
outdir <- snakemake@params[['outdir']]
report_fp <- snakemake@output[['report']]

## Parse
blastn <- do.call(rbind, 
                  lapply(1:length(summary_files), 
                         function(x) read_delim(summary_files[x], delim="\t"))) %>% 
  filter(length >= 2000) %>%
  filter_(paste("!is.na(", db, ")", sep="")) %>%
  mutate(taxaID = accessionToTaxa(get(db), accessionTaxasql)) 

## Filter by taxaName, this is slower, yet RefSeq doesn't necessarily update the new assembly 
blastn <- cbind(blastn, getTaxonomy(blastn$taxaID,taxaNodes,taxaNames))

blastn_taxa <- filter(blastn, grepl(taxaName, species, perl=T, ignore.case = T))

## NOW all we need is to write the contig names of the filtered contigs into separate file

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
