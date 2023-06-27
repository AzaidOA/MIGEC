library(stringr)

args = commandArgs(trailingOnly=TRUE)

migec.work.dir <- args[1]

sample.name <- basename(getwd())
migec.checkout.dir <- paste0(migec.work.dir,'/',sample.name,"/checkout")
migec.histogram.dir <- paste0(migec.work.dir,'/',sample.name,"/histogram")

if(!dir.exists(migec.histogram.dir)) dir.create(migec.histogram.dir)

line <- paste(migec.checkout.dir, migec.histogram.dir, sep = " ")

cat(line)
