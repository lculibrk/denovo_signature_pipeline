' tsv_to_maf.R

Usage: tsv_to_maf.R -i INPUT -c COLUMNS -o OUTPUT [ --allow-duplicates ]

Options:
    -i --input INPUT        Path to input TSV
    -c --columns COLUMNS    Comma-separated list of sample ID, chromosome, position, ref, and alt columns in that order
    -o --output OUTPUT      Path to output
    --allow-duplicates      Include this flag if you wish to allow lines with duplicate loci
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

snv <- read_tsv(args[['input']], col_types = cols(.default='c'))
columns = strsplit(args[['columns']], ',')[[1]]

snv <- snv[, columns]
colnames(snv) <- c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Allele')

snv <- snv %>%
    group_by(Tumor_Sample_Barcode, Chromosome, Start_Position) %>%
    summarise(
        Reference_Allele = as.character(names(sort(table(Reference_Allele), decreasing=TRUE)[1])),
        Allele = as.character(names(sort(table(Allele), decreasing=TRUE)[1]))
    )

write_tsv(snv, args[['output']])
