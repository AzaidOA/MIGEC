# MIGEC
In-house developed pipeline for bulk TCR data preprocessing based on the Molecular Identifier-Guided Error Correction (MIGEC) tool.<br/>
MIGEC was [published](https://www.nature.com/articles/nmeth.2960) and maintained by [Dr. Mikhail Shugay](https://github.com/mikessh). Our in-house pipeline was intended to help process the data from hundreds or even couples of thousands of bulk TCR libraries in an efficient manner.

## Introduction

The MIGEC tool is designed to process molecular barcoding data, correct errors, and generate a consensus sequence for each unique molecular identifier. It is widely used in the field of single-cell sequencing and immune repertoire analysis. MIGEC helps improve the accuracy of downstream analyses by addressing the challenges posed by errors introduced during sequencing and PCR amplification.

For more information on MIGEC, refer to the following resources:
- [MIGEC Paper](https://www.nature.com/articles/nmeth.2960)
- [MIGEC GitHub Repository](https://github.com/mikessh/migec)

## Requirements

To run the MIGEC pipeline, you will need the following:

- R (version 3.6.2 or higher)
- R packages:
  - data.table
  - stringr
  - plyr
  - reshape2

Additionally, you will need to install NCBI-BLAST+. Please follow these steps to install it:

1. Visit the [NCBI-BLAST+ website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to download the appropriate package for your operating system.

2. Follow the installation instructions provided on the website to install NCBI-BLAST+. Make sure to install it in a location that is accessible from the command line.

## Installation

To set up the MIGEC pipeline, follow these steps:

1. Clone or download this repository to your local machine.

2. Navigate to the `R` directory inside the repository.

3. Make sure you have the required R packages installed: `data.table`, `stringr`, `plyr`, `reshape2`. If you don't have them installed, you can install them by running the following commands in your R console:

 ```R
 install.packages("data.table")
 install.packages("stringr")
 install.packages("plyr")
 install.packages("reshape2")
 ```

 4. Run the run_migec.R script with the provided command-line arguments. Adjust the arguments according to your setup and requirements.

--ReportsDir: Absolute path to the directory to save reports.
--ProjectID: Name of the project to run.
--bcl2fastqSheet: Absolute path to the sample sheet used to run bcl2fastq.
--bcl2fastqReports: Absolute path to the directory with the reports of bcl2fastq.
--BarcodesSheet: Absolute path to the sheet with the master and slave barcodes.
--TimeLimit: Maximum minutes to allow the pipeline to aggregate all the fastq files.
--TimeLimitMIGEC: Maximum minutes to allow the MIGEC tool to run.
--Species: Default species to be used for all samples. Options: 'HomoSapiens', 'MusMusculus'.
--FileType: File type to be processed for all the samples. Options: 'paired', 'overlapped', 'single'.
--Mask: Mask which specifies for which reads in paired-end data to perform the CDR3 extraction. R1='1:0', R2='0:1', Both='1:1'.

Rscript run_migec.R --ReportsDir /path/to/reports --ProjectID my_project --bcl2fastqSheet /path/to/sample_sheet.csv --bcl2fastqReports /path/to/bcl2fastq_reports --BarcodesSheet /path/to/barcodes_sheet.csv --TimeLimit 120 --TimeLimitMIGEC 300 --Species HomoSapiens --FileType paired --Mask 1:1

Feel free to customize any sections according to your specific needs. Let me know if you need any further assistance!
