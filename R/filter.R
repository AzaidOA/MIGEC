library(stringr)

args = commandArgs(trailingOnly=TRUE)

migec.work.dir <- args[1]

sample.name <- basename(getwd())
migec.cdrblast.dir <- paste0(migec.work.dir,'/', sample.name, "/cdrblast")
migec.filter.dir <- paste0(migec.work.dir, '/', sample.name, "/filter")

if(!dir.exists(migec.filter.dir)) dir.create(migec.filter.dir)
line <- paste(migec.cdrblast.dir, migec.filter.dir)

cat(line)
