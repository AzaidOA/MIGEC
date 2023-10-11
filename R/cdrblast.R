library(stringr)

args = commandArgs(trailingOnly=TRUE)

migec.work.dir <- args[1]

sample.name <- basename(getwd())
migec.checkout.dir <- paste0(migec.work.dir, '/', sample.name, "/checkout")
migec.assembly.dir <- paste0(migec.work.dir, '/', sample.name, "/assembly")
migec.cdrblast.dir <- paste0(migec.work.dir, '/', sample.name, "/cdrblast")

if(!dir.exists(migec.cdrblast.dir)) dir.create(migec.cdrblast.dir)

metadata.file <- list.files(path = getwd(), pattern = "metadata", full.names = T)

line <- paste("--sample-metadata", metadata.file, migec.checkout.dir, migec.assembly.dir, migec.cdrblast.dir)

cat(line)
