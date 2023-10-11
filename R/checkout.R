library(stringr)

args = commandArgs(trailingOnly=TRUE)

migec.work.dir <- args[1]
migec.data.dir <- getwd()
migec.checkout.dir <- paste0(migec.work.dir, '/',
                               basename(migec.data.dir))
if(!dir.exists(migec.checkout.dir)) dir.create(migec.checkout.dir)
migec.checkout.dir <- paste0(migec.checkout.dir,"/checkout")
if(!dir.exists(migec.checkout.dir)) dir.create(migec.checkout.dir)

files <- list.files(migec.data.dir)
line <- paste(files[str_detect(files, "barcodes")],
              files[str_detect(files, "R1")], files[str_detect(files, "R2")],
              migec.checkout.dir)
cat(line)
