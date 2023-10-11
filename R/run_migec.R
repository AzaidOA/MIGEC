# Main: MIGEC pipeline
# By: Azaid Ordaz
# Version: 1.1.0

cat('############    --------   Bulk TCR-seq   -------    ############\n')
cat('############    ------------   MIGEC   -----------    ############\n')
### ------------------------------ Libraries ------------------------------ ###
cat('### ------------------------------ Libraries ------------------------------ ###\n')
cat('Importing libraries...\n')
library(optparse)
library(data.table)
library(stringr)
library(plyr)
library(reshape2)
### ------------------------------ Arguments ------------------------------ ###
cat('### ------------------------------ Arguments ------------------------------ ###\n')
option.list <- list(
  make_option(opt_str='--ReportsDir', type='character', default=NULL, dest='migec.work.dir', help='Absolute path to directory to save reports.'),
  make_option(opt_str='--ProjectID', type='character', default=NULL, dest='project.id', help='Name of the project to run.'),
  make_option(opt_str='--bcl2fastqSheet', type='character', default=NULL, dest='bcl2fastq.sample.sheet.file', help='Absolute path to sample sheet used to run bcl2fastq.'),
  make_option(opt_str='--bcl2fastqReports', type='character', default=NULL, dest='bcl2fastq.reports.dir', help='Absolute path to the directory with the reports of bcl2fastq.'),
  make_option(opt_str='--BarcodesSheet', type='character', default=NULL, dest='barcodes.file', help='Absolute path to sheet with the master and slave barcodes.'),
  make_option(opt_str='--SampleSheet', type='character', default=NULL, dest='migec.sample.sheet.file', help='Absolute path to the sheet with the data related to each sample. Needed columns: \'sample.id\', \'library.id\', \'donor.id.tag\', \'chain.tag\', \'master.id\', \'slave.id\'.'),
  make_option(opt_str='--TimeLimit', type='character', default='06:00:00', dest='time.limit', help='Maximum time in format \'06:00:00\' \'hours:minutes:seconds\' to allow MIGEC to run.'),
  make_option(opt_str='--Species', type='character', default='HomoSapiens', dest='def.specie', help='Default species to be used for all samples, \'HomoSapiens\' and \'MusMusculus\' available.'),
  make_option(opt_str='--FileType', type='character', default='paired', dest='def.file.type', help='File type to be processed for all the samples, \'paired\', \'overlapped\' and \'single\' available.'),
  make_option(opt_str='--Mask', type='character', default='1:1', dest='def.mask', help='Mask which specifies for which reads in paired-end data to perform the CDR3 extraction. R1=\'1:0\', R2=\'0:1\', Both=\'1:1\', and Only Overlapped reads=\'0:0\''),
  make_option(opt_str='--Quality', type='character', default='25,30', dest='def.quality', help='Quality threshold pair, in comma-separated format (i.e. \'25,30\'), default for all samples. First threshold in pair is used for raw sequence quality (sequencing quality phred) and the second one is used for assembled sequence quality (CQS score, the fraction of reads in MIG that contain dominant letter at a given position).'),
  make_option(opt_str='--MIGsize', type='integer', default='2', dest='MIG.size', help='Minimum number of reads that each UMI must have.')
  
)


if(interactive()){

  project.id <- 'BulkTCR016'
  migec.work.dir <- '/mnt/BioAdHoc/Groups/vd-vijay/moarias/trial/MIGEC'
  bcl2fastq.sample.sheet.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/sequencing_data/04-17-2023/mkfastq/data/experiment_layout/NV103_sample_sheet_for_mkfastq.csv'
  bcl2fastq.reports.dir <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/sequencing_data/04-17-2023/mkfastq/NV103/outs/fastq_path/BulkTCR016'
  barcodes.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC/metadata/BulkTCR016_barcodes_sheet_migec.csv'
  migec.sample.sheet.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC/metadata/BulkTCR016_sample_sheet_migec.csv'
  time.limit <- '06:00:00'
  def.specie <- "HomoSapiens"
  def.file.type <- "paired"
  def.mask <- "1:1"
  def.quality <- "25,30"
  MIG.size <- 2


}

opt.parser <- OptionParser(option_list=option.list)
opt <- parse_args(opt.parser)

project.id <- opt$project.id
migec.work.dir <- opt$migec.work.dir
bcl2fastq.sample.sheet.file <- opt$bcl2fastq.sample.sheet.file
bcl2fastq.reports.dir <- opt$bcl2fastq.reports.dir
barcodes.file <- opt$barcodes.file
migec.sample.sheet.file <- opt$migec.sample.sheet.file
time.limit <- opt$time.limit
def.specie <- opt$def.specie
def.file.type <- opt$def.file.type
def.mask <- opt$def.mask
def.quality <- opt$def.quality
MIG.size <- opt$MIG.size


# ---> Sample sheet used to run bcl2fastq.
if(!file.exists(bcl2fastq.sample.sheet.file))
  stop('bcl2fastq sample sheet was not found.\n')

# ---> Directory with the results from bcl2fastq.
if(!dir.exists(bcl2fastq.reports.dir))
  stop('The path to the directory with the bcl2fastq results was not found.\n')


# ---> Barcodes sheet.
if(is.null(barcodes.file))
  barcodes.file <- paste0(migec.work.dir,'/metadata/',project.id,'_barcodes_sheet_migec.csv')
if(!file.exists(barcodes.file))
  stop('Barcodes sheet for MIGEC was not found.\n')

# ---> Metadata sheet.
if(is.null(migec.sample.sheet.file))
  migec.sample.sheet.file <- paste0(migec.work.dir,'/metadata/',project.id,'_sample_sheet_migec.csv')
if(!file.exists(migec.sample.sheet.file))
  stop('Sample sheet for MIGEC was not found.\n')

# ---> Directory where the data will be saved.
migec.data.dir <- paste0(migec.work.dir, '/data')
if(!dir.exists(migec.data.dir)) dir.create(migec.data.dir)
migec.data.dir <- paste0(migec.data.dir, '/', project.id)
if(!dir.exists(migec.data.dir)) dir.create(migec.data.dir)

# ---> Directory where the jobs scripts will be saved.
migec.scripts.dir <- paste0(migec.work.dir,'/jobs_scripts')
if(!dir.exists(migec.scripts.dir)) dir.create(migec.scripts.dir)
migec.scripts.dir <- paste0(migec.scripts.dir,'/', project.id)
if(!dir.exists(migec.scripts.dir)) dir.create(migec.scripts.dir)

# ---> Directory where MIGEC will run.
migec.outs.dir <- paste0(migec.work.dir,"/outs")
if(!dir.exists(migec.outs.dir)) dir.create(migec.outs.dir)
migec.outs.dir <- paste0(migec.outs.dir,'/',project.id)
if(!dir.exists(migec.outs.dir)) dir.create(migec.outs.dir)

### ------------------- Joining the fastqs of each Read ------------------- ###
cat('### -------------------- Agreggating fastq files -------------------- ###\n')
bcl2fastq.sample.sheet <- fread(bcl2fastq.sample.sheet.file,skip=1)
# Getting library names.
libraries <- bcl2fastq.sample.sheet[Sample_Project == project.id,unique(Original_Sample_ID)]

# Creating one directory per library.
data.dirs <- paste0(migec.data.dir,'/',libraries)
invisible(sapply(X=data.dirs,FUN=dir.create))
cat('Data directories were created...\n\n')

# Getting fastq files.
fastq.files <- list.files(path=bcl2fastq.reports.dir, recursive=T, full.names=T)

aggregate.fastq.file <- paste0(migec.scripts.dir,'/',project.id,'_aggregate_fastqs.sh')
aggregate.fastq.flag.file <- paste0(migec.scripts.dir,'/aggregate_fastq_files_done.txt')
aggregate.fastq <- c()

for(tmp.library in libraries){
  r1 <- paste0(fastq.files[str_detect(fastq.files,paste0(tmp.library,'_.*R1.+'))],collapse=' ')
  r2 <- paste0(fastq.files[str_detect(fastq.files,paste0(tmp.library,'_.*R2.+'))],collapse=' ')

  aggregate.fastq <- c(aggregate.fastq,
    paste0('cat ',r1,' > ',migec.data.dir,'/',tmp.library,'/',tmp.library,'_R1.fastq.gz'),
    paste0('cat ',r2,' > ',migec.data.dir,'/',tmp.library,'/',tmp.library,'_R2.fastq.gz')
  )
}

write(x=aggregate.fastq, file=aggregate.fastq.file)
system(command=paste0('bash ', aggregate.fastq.file))
cat('FASTQ file were agreggated succesfully...\n')

### ----------------------- Creating barcodes files ----------------------- ###
cat('### ----------------------- Creating barcodes files ----------------------- ###\n')
migec.sample.sheet <- fread(file=migec.sample.sheet.file)
barcodes.data <- fread(barcodes.file)
barcodes.files.list <- migec.sample.sheet[,.(library.id,sample.id,master.id,slave.id)]
barcodes.files.list[ barcodes.data[type=='master'],
  on=.(master.id=barcode.id), master.sequence:=sequence]
barcodes.files.list[ barcodes.data[type=='slave'],
  on=.(slave.id=barcode.id), slave.sequence:=sequence]

barcodes.files.list[,':='(
  master.id=NULL, slave.id=NULL,
  master.sequence=paste0(master.sequence,'tNNNNtNNNNtNNNNtct'),
  slave.sequence=paste0("NNNN",slave.sequence))]

colnames(barcodes.files.list) <- c("library.id","SAMPLE ID", "Master barcode sequence", "Slave barcode sequence")
barcodes.files.list <- split(x=barcodes.files.list, by="library.id", keep.by=F)

migec.libraries <- unique(migec.sample.sheet$library.id)
for(tmp.lib in migec.libraries){
  tmp.file <- paste0(migec.data.dir,'/',tmp.lib,"/barcodes_",tmp.lib,".txt")
  fwrite(x=barcodes.files.list[[tmp.lib]], file=tmp.file, sep="\t")
}

cat('Barcodes files were created succesfully...\n')

### ----------------------- Creating metadata files ----------------------- ###
cat('### ----------------------- Creating metadata files ----------------------- ###\n')
migec.sample.sheet <- fread(file=migec.sample.sheet.file)

metadata.files.list <- migec.sample.sheet[,.(
  library.id,sample.id,def.specie,chain.tag,def.file.type,def.mask,def.quality)]
colnames(metadata.files.list) <- c("library.id","Sample ID","Species","Gene","File types","Mask","Quality threshold pair")
metadata.files.list <- split(x=metadata.files.list, by="library.id", keep.by=F)

migec.libraries <- unique(migec.sample.sheet$library.id)
for(tmp.lib in migec.libraries){
  tmp.file <- paste0(migec.data.dir,'/',tmp.lib,"/metadata_",tmp.lib,".txt")
  fwrite(x=metadata.files.list[[tmp.lib]], file=tmp.file, sep="\t", col.names=F)
}
cat('Metadata files were created succesfully...\n')


### -------------------------- Creating job files -------------------------- ###
cat('### -------------------------- Creating job files -------------------------- ###')
migec.libraries <- list.files(path=migec.data.dir, recursive=F)
SBATCH.dir <- paste0(migec.scripts.dir,"/SBATCH")
if(!dir.exists(SBATCH.dir)) dir.create(SBATCH.dir)

run.migec <- c()
for(tmp.lib in migec.libraries){
  SBATCH <- paste0("#!/bin/sh
#SBATCH --job-name=MIGEC-",tmp.lib,"
#SBATCH --time=",time.limit,"
#SBATCH --output=", paste0(migec.data.dir,"/",tmp.lib,"/MIGEC.out\n"),
"#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
WORK_DIR=\"",migec.outs.dir,"\"
SAMPLE_DIR=\"",paste0(migec.data.dir,"/",tmp.lib), "\"
MIGEC=\"java -jar /home/moarias/bin/MIGEC/migec-1.2.9.jar\"

cd ${SAMPLE_DIR}
#Running checkout
checkout_line=$(Rscript /home/moarias/scripts/migec_processing/checkout.R \\
${WORK_DIR})
${MIGEC} Checkout -cute ${checkout_line}

#Running histogram
histogram_line=$(Rscript /home/moarias/scripts/migec_processing/histogram.R \\
${WORK_DIR})
${MIGEC} Histogram ${histogram_line}

#Running assemble
assembly_line=$(Rscript /home/moarias/scripts/migec_processing/assembly.R \\
${WORK_DIR})
${MIGEC} AssembleBatch --force-overseq ",MIG.size," --force-collision-filter -c ${assembly_line}

#Running cdrblastbatch
cdrblast_line=$(Rscript /home/moarias/scripts/migec_processing/cdrblast.R \\
${WORK_DIR})
${MIGEC} CdrBlastBatch --all-alleles ${cdrblast_line}

#Running filterning
filter_line=$(Rscript /home/moarias/scripts/migec_processing/filter.R \\
${WORK_DIR})
${MIGEC} FilterCdrBlastResultsBatch -s ${filter_line}")
  tmp.file <- paste0(SBATCH.dir,"/SBATCH_",tmp.lib,".sh")
  run.migec <- c(run.migec,tmp.file)
  write(x=SBATCH, file=tmp.file)
}

run.migec <- paste("sbatch", run.migec)
write(x=run.migec,file=paste0(migec.scripts.dir,"/run_migec.sh"))
system(command=paste0('bash ', migec.scripts.dir,"/run_migec.sh"))

cmd.files <- paste0(migec.outs.dir,'/',migec.libraries,'/filter/cdrblastfilter.cmd.txt')
names(cmd.files) <- migec.libraries
all.out.files <- paste0(migec.data.dir,'/',migec.libraries,'/MIGEC.out')
migec.done <- all(file.exists(cmd.files))
time.limit <- str_extract_all(time.limit, "\\d+")[[1]]
time.limit <- (as.numeric(time.limit[1])*60) + as.numeric(time.limit[2])

time.count <- 0
errors <- list()
while(!migec.done){
  # Checking if any errors have arisen...
  out.files <- list.files(path=migec.data.dir, full.names=T, recursive=T, pattern='MIGEC.out')
  names(out.files) <- basename(dirname(out.files))
  # Removing MIGEC.out files that have already been detected as an error...
  out.files <- out.files[!names(out.files) %chin% names(errors)]
  new.errors <- lapply(X=out.files,FUN=function(tmp.file){
    grep.command <- paste("grep", 'ERROR', shQuote(tmp.file), sep=" ")
    tmp.errors <- system(grep.command, intern=TRUE)
    if(length(tmp.errors) == 0) return(NA)
    return(tmp.errors)
  })
  new.errors <- new.errors[!is.na(new.errors)]
  
  if(length(new.errors) > 0) cmd.files <- cmd.files[!names(cmd.files) %chin% names(new.errors)]

  # Error detection: Division undefined
  low.umis <- sapply(new.errors,FUN=function(tmp.error){
    any(tmp.error %chin% '[ERROR] Division undefined, see _migec_error.log for details')
  })
  if(length(low.umis) > 0){
    new.errors <- new.errors[!low.umis]
    low.umis <- names(low.umis[low.umis])
    for(tmp.lib in low.umis){
      # Detecting which samples were affected
      assemble.file <- paste0(migec.outs.dir,'/',tmp.lib, '/assembly/assemble.log.txt')
      assemble <- fread(assemble.file, select=c('#SAMPLE_ID','MIGS_GOOD_TOTAL'))
      samples.to.remove <- as.character(unlist(assemble[MIGS_GOOD_TOTAL == 0, '#SAMPLE_ID']))

      # Removing samples from metadata and barcode files
      barcode.file <- paste0(migec.data.dir,'/',tmp.lib,'/barcodes_',tmp.lib,'.txt')
      barcodes <- fread(barcode.file, colClasses='character')
      barcodes <- barcodes[!`SAMPLE ID` %chin% samples.to.remove]
      fwrite(x=barcodes, file=barcode.file, sep="\t")

      metadata.file <- paste0(migec.data.dir,'/',tmp.lib,'/metadata_',tmp.lib,'.txt')
      metadata <- fread(metadata.file, colClasses='character')
      metadata <- metadata[!V1 %chin% samples.to.remove]
      fwrite(x=metadata, file=metadata.file, col.names=F, sep="\t")

      # Creatig log file
      samples.to.remove <- paste(samples.to.remove, collapse=' ')
      log.message <- paste0('The following samples were eliminated because they did not reach the requested MIG size ', MIG.size,': ',samples.to.remove)
      log.file <- paste0(migec.data.dir,'/',tmp.lib,'/removed_samples.log.txt')
      write(x=log.message, file=log.file)

      # Remiving MIGEC.out and _migec_error.log files
      files.to.remove <- c(
        paste0(migec.data.dir,'/',tmp.lib,'/MIGEC.out'), 
        paste0(migec.data.dir,'/',tmp.lib,'/_migec_error.log')
      )
      invisible(file.remove(files.to.remove))

      # Keeping histogram log.
      histogram.dir <- paste0(migec.outs.dir,'/',tmp.lib, '/histogram')
      histogram.dir.log <- paste0(migec.outs.dir,'/',tmp.lib, '/log_histogram')
      invisible(file.rename(histogram.dir, histogram.dir.log))

      # Removing all the other directories.
      dirs.to.remove <- list.dirs(path=paste0(migec.outs.dir,'/',tmp.lib), recursive=F)
      dirs.to.remove <- dirs.to.remove[basename(dirs.to.remove) != 'log_histogram']
      invisible(unlink(dirs.to.remove, recursive=T))

      # Re-runing job.
      sbatch.line <- paste0('sbatch ', SBATCH.dir, '/SBATCH_', tmp.lib,'.sh')
      system(command=sbatch.line)

      cat(paste0('The samples from the ',tmp.lib,' library that did not reach the defined MIG size ',MIG.size,' have been eliminated and the job has been run again.\n'))
    }
    time.count <- 0
  }
  
  errors <- append(errors, new.errors)
  time.start <- all(file.exists(all.out.files))
  if(time.start){
    if(time.count > time.limit){
    stop('MIGEC did not run in the allocated time. Check for any incompleted library.\n')
    }
  time.count <- time.count + 1
  }
  Sys.sleep(time=60)
  migec.done <- all(file.exists(cmd.files))
}
if(length(errors) > 0){
  project.error.file <- paste0(migec.data.dir,'/project_error.log.txt')
  errors <- lapply(errors, paste, collapse='\n')
  log.conn <- file(project.error.file, "w")
  for(i in 1:length(errors)){
    cat(paste0(names(errors)[i],'\n'), file=log.conn)
    cat(paste0(errors[[i]],'\n\n'), file=log.conn)
  }
  close(log.conn)
  cat('The following libraries had errors while running, check the specific errors they got in project_error.log.txt file:\n')
  cat(paste(names(errors),collapse='\n'))
  cat('\nBut the other libraries have run successfully!\n')
} else {
  cat('MIGEC ran succesfully!\n')
}
 

### -------------------------- Creating summaries -------------------------- ###
cat('### -------------------------- Creating summaries -------------------------- ###\n')
library.dirs <- list.dirs(migec.outs.dir, full.names = T, recursive = F)
#pool.dirs <- pool.dirs[basename(pool.dirs) %like% 'A']

filters <- c("raw","asm")
## Getting the sample files per pool
for (f in filters){
  for (tmp.dir in library.dirs){
    cdrblast.dir <- paste0(tmp.dir, "/cdrblast")
    cdrblast.files <- list.files(cdrblast.dir, full.names = T, pattern = paste0(f,".cdrblast"))
    ## Geeting the information from the clonotyope tables of each sample
    library.summary <- data.table()
    for (sample.file in cdrblast.files){
      clonotypes.table <- fread(sample.file,
                            select = c("Count", "CDR3 amino acid sequence", "V segments", "J segments", "Good events", "Good reads"),
                            col.names = c("count", "cdr3.sequence", "v.gene", "j.gene", "total.umis", "total.reads"))
      ## Getting sumarry of the clonotype table
      sample.id <- str_remove(basename(sample.file), paste0(".",f,".+$"))
      total.clonotypes <- nrow(clonotypes.table)
      total.reads <- clonotypes.table[,sum(total.reads)]
      total.umis <- clonotypes.table[,sum(total.umis)]
      mean.clone.size <- clonotypes.table[,mean(count)]
      library.summary <- rbind(library.summary,
                            data.table(sample.id, sample.file, total.reads, total.umis, total.clonotypes, mean.clone.size))
    }
    fwrite(x = library.summary, file = paste0(tmp.dir, "/cdrblast/library.",f,".summary.txt"), sep = "\t")
  }
}

## Concatenating asm and raw if exist
cdrblast.dirs <- paste0(library.dirs, "/cdrblast")
for (tmp.dir in cdrblast.dirs){
  cdrblast.files <- list.files(tmp.dir, full.names = T, pattern = "^library.(raw|asm)")
  cdrblast.summary <- rbind(
    fread(cdrblast.files[1])[, data.type := str_extract(basename(cdrblast.files[1]), "raw|asm")],
    fread(cdrblast.files[2])[,data.type := str_extract(basename(cdrblast.files[2]), "raw|asm")]
  )
  fwrite(x = cdrblast.summary, file = paste0(tmp.dir, "/library.cdrblast.summary.txt"), sep = "\t")
}

## Filter

for (tmp.dir in library.dirs){
  filter.dir <- paste0(tmp.dir, "/filter")
  filter.files <- list.files(filter.dir, full.names = T, pattern = paste0(".filtered.cdrblast.txt"))
  ## Geeting the information from the clonotyope tables of each sample
  library.summary <- data.table()
  for (sample.file in filter.files){
    clonotypes.table <- fread(sample.file,
                              select = c("Count", "CDR3 amino acid sequence", "V segments", "J segments", "Good events", "Good reads"),
                              col.names = c("count", "cdr3.sequence", "v.gene", "j.gene", "total.umis", "total.reads"))
    ## Getting sumary of the clonotype table
    sample.id <- str_remove(basename(sample.file), paste0(".filtered.cdrblast.txt$"))
    total.clonotypes <- nrow(clonotypes.table)
    total.reads <- clonotypes.table[,sum(total.reads)]
    total.umis <- clonotypes.table[,sum(total.umis)]
    mean.clone.size <- clonotypes.table[,mean(count)]
    library.summary <- rbind(library.summary,
                          data.table(sample.id, sample.file, total.reads, total.umis, total.clonotypes, mean.clone.size))
  }
  fwrite(x = library.summary, file = paste0(tmp.dir, "/filter/library.filter.summary.txt"), sep = "\t")
}

# ---> Creating QCs

summary.dir <- paste0(migec.work.dir,"/summaries")
if(!dir.exists(summary.dir)) dir.create(summary.dir)
summary.dir <- paste0(summary.dir,'/',project.id)
if(!dir.exists(summary.dir)) dir.create(summary.dir)
summary.file <- paste0(summary.dir,"/summary_per_sample.csv")

migec.sample.sheet <- fread(file=migec.sample.sheet.file, drop=c("master.id","slave.id"))

migec.sample.sheet[,':='(unfiltered.reads = 0, unfiltered.umis = 0,
                   first.reads = 0, first.umis = 0,
                   first.clonotypes = 0, first.mean.clone.size = 0,
                   second.reads = 0, second.umis = 0,
                   second.clonotypes = 0, second.mean.clone.size = 0)]

out.directories <- list.dirs(migec.outs.dir, full.names = T, recursive = F)

for (tmp.dir in out.directories) {
 ## cdrblast.summary: The metrics of the first filter step that I created from the clonotype tables
 ## cdrblast.log: Metrics of the first filter step directly created by MIGEC
 ## filter.summary; Metrics of the second filter step created by me, but are the same that MIGEC creates
 filter.summary.file <- paste0(tmp.dir, "/filter/library.filter.summary.txt")
 cdrblast.summary.file <- paste0(tmp.dir,"/cdrblast/library.cdrblast.summary.txt")
 cdrblast.log.file <- paste0(tmp.dir,"/cdrblast/cdrblast.log.txt")
 filter.summary.table <- fread(filter.summary.file, drop = "sample.file")
 cdrblast.summary.table <- fread(cdrblast.summary.file, drop = "sample.file")[data.type == "asm"]
 cdrblast.log.table <- fread(cdrblast.log.file, select = c("#SAMPLE_ID", "DATA_TYPE", "READS_GOOD", "EVENTS_GOOD", "READS_TOTAL", "EVENTS_TOTAL"),
                             col.names = c("sample.id", "data.type", "cdrblast.total.reads", "cdrblast.total.umis", "no.cdrblast.total.reads", "no.cdrblast.total.umis"))
 library.d <- basename(tmp.dir)
 for (id.r in filter.summary.table$sample.id){
   migec.sample.sheet[sample.id == id.r & library.id == library.d,
   ':='(unfiltered.reads = cdrblast.log.table[sample.id == id.r & data.type == "raw", no.cdrblast.total.reads],
   unfiltered.umis = cdrblast.log.table[sample.id == id.r & data.type == "raw", no.cdrblast.total.umis],
   first.reads = cdrblast.summary.table[sample.id == id.r, total.reads],
   first.umis = cdrblast.summary.table[sample.id == id.r, total.umis],
   first.clonotypes = cdrblast.summary.table[sample.id == id.r, total.clonotypes],
   first.mean.clone.size = cdrblast.summary.table[sample.id == id.r, mean.clone.size],
   second.reads = filter.summary.table[sample.id == id.r, total.reads],
   second.umis = filter.summary.table[sample.id == id.r, total.umis],
   second.clonotypes = filter.summary.table[sample.id == id.r, total.clonotypes],
   second.mean.clone.size = filter.summary.table[sample.id == id.r, mean.clone.size])]
  }
}

fwrite(x = migec.sample.sheet, file = summary.file)
