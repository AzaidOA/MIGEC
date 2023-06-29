# Main: MIGEC pipeline
# By: Azaid Ordaz
# Version: 1.0.0

cat('############    --------   Bulk TCR-seq   -------    ############\n')
cat('############    ------------   MIGEC   -----------    ############\n')
### ------------------------------ Libraries ------------------------------ ###
cat('### ------------------------------ Libraries ------------------------------ ###\n')
cat('Importing libraries...\n')
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
  make_option(opt_str='--TimeLimit', type='integer', default=300, dest='time.limit', help='Maximum minutes to allow MIGEC to run.'),
  # make_option(opt_str='--TimeLimitMIGEC', type='integer', default=300, dest='migec.time.limit', help='Time limit (in minutes) for MIGEC to run.'),
  make_option(opt_str='--Species', type='character', default='HomoSapiens', dest='def.specie', help='Default species to be used for all samples, \'HomoSapiens\' and \'MusMusculus\' available.'),
  make_option(opt_str='--FileType', type='character', default='paired', dest='def.file.type', help='File type to be processed for all the samples, \'paired\', \'overlapped\' and \'single\' available.'),
  make_option(opt_str='--Mask', type='character', default='1:1', dest='def.mask', help='Mask which specifies for which reads in paired-end data to perform the CDR3 extraction. R1=\'1:0\', R2=\'0:1\', Both=\'1:1\', and Only Overlapped reads=\'0:0\''),
  make_option(opt_str='--Quality', type='character', default='25,30', dest='def.quality', help='Quality threshold pair, in comma-separated format (i.e. \'25,30\'), default for all samples. First threshold in pair is used for raw sequence quality (sequencing quality phred) and the second one is used for assembled sequence quality (CQS score, the fraction of reads in MIG that contain dominant letter at a given position).')
)


if(interactive()){

  # BulkTCR016
  # seq.date <- '04-17-2023'
  # run.id <- 'NV103'
  project.id <- 'BulkTCR016'
  migec.work.dir <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/trial/MIGEC'
  bcl2fastq.sample.sheet.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/sequencing_data/04-17-2023/mkfastq/data/experiment_layout/NV103_sample_sheet_for_mkfastq.csv'
  bcl2fastq.reports.dir <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/sequencing_data/04-17-2023/mkfastq/NV103/outs/fastq_path/BulkTCR016'
  barcodes.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC/metadata/BulkTCR016_barcodes_sheet_migec.csv'
  migec.sample.sheet.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC/metadata/BulkTCR016_sample_sheet_migec.csv'
  time.limit <- 300
  def.specie <- "HomoSapiens"
  def.file.type <- "paired"
  def.mask <- "1:1"
  def.quality <- "25,30"

  # Bulk_TCR008
  seq.date <- '04-05-2022'
  run.id <- 'NV078'
  project.id <- 'BulkTCR_08-10-11-12'
  migec.work.dir <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC'
  time.limit <- 120
  master.extension <- "cagtggtatcaacgcagagtNNNNtNNNNtNNNNtct"
  trb.slave.extension <- "acacsttkttcaggtcctc"
  tra.slave.extension <- "gggtcagggttctggatat"

}

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
# It has to exist project directory, and inside all its respective fastq files
# aggregate.fastq <- c(paste0('cat ',bcl2fastq.outs.dir,'/',libraries,'*R1* > ',migec.data.dir,'/',project.id,'/',libraries,'/',libraries,'_R1.fastq.gz'),
#   paste0('cat ',bcl2fastq.outs.dir,'/',libraries,'*R2* > ',migec.data.dir,'/',project.id,'/',libraries,'/',libraries,'_R2.fastq.gz'))
# RIN
# aggregate.fastq <- c(aggregate.fastq,
#   paste0('echo \'FASTQ files are ready. This file has to be automatically removed\' > ',
#     aggregate.fastq.flag.file))
# RIN
# batch <- project.id
write(x=aggregate.fastq, file=aggregate.fastq.file)
system(command=paste0('bash ', aggregate.fastq.file))
cat('FASTQ file were agreggated succesfully...\n')
# RIN
# time.count <- 0
# while(!file.exists(aggregate.fastq.flag.file)){
#   time.count <- time.count + 1
#   if(time.count > time.limit){
#     stop('The FASTQ files were not aggregated in the allocated time.')
#   }
#   Sys.sleep(time=60)
# }
# system(command=paste0('rm ', aggregate.fastq.flag.file))
### ----------------------- Creating barcodes files ----------------------- ###
cat('### ----------------------- Creating barcodes files ----------------------- ###\n')
migec.sample.sheet <- fread(file=migec.sample.sheet.file)
#barcodes.files.list <- migec.sample.sheet[,.(sample.id,master.sequence,slave.sequence), by="pool"]
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


### -------------------------- Creating job files -------------------------- ###
cat('### -------------------------- Creating job files -------------------------- ###')
# project.data.dir <- paste0(migec.data.dir,'/',project.id)
migec.libraries <- list.files(path=migec.data.dir, recursive=F)
#migec.libraries <- migec.libraries[migec.libraries %like% "A"]
SBATCH.dir <- paste0(migec.scripts.dir,"/SBATCH")
if(!dir.exists(SBATCH.dir)) dir.create(SBATCH.dir)

run.migec <- c()
for(tmp.lib in migec.libraries){
  SBATCH <- paste0("#!/bin/sh
#SBATCH --job-name=MIGEC-",tmp.lib,"
#SBATCH --time=",time.limit/60,":00:00
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
${MIGEC} AssembleBatch --force-overseq 2 --force-collision-filter -c ${assembly_line}

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
migec.done <- all(file.exists(cmd.files))
time.count <- 0
while(!migec.done){
  time.count <- time.count + 1
  migec.done <- all(file.exists(cmd.files))
  if(time.count > time.limit){
    stop('MIGEC did not run in the allocated time. Check for any incompleted library.\n')
  }
  Sys.sleep(time=60)
}
cat('MIGEC ran succesfully!\n')

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
