library(stringr)

args = commandArgs(trailingOnly=TRUE)

migec.work.dir <- args[1]


sample.name <- basename(getwd())
migec.checkout.dir <- paste0(migec.work.dir,'/',sample.name,"/checkout")
migec.histogram.dir <- paste0(migec.work.dir,'/',sample.name,"/histogram")
migec.assembly.dir <- paste0(migec.work.dir,'/',sample.name,"/assembly")


if(!dir.exists(migec.assembly.dir)) dir.create(migec.assembly.dir)

line <- paste(migec.checkout.dir, migec.histogram.dir, migec.assembly.dir)

cat(line)
