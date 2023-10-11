#!/bin/sh
#SBATCH --job-name=MIGEC-BulkTCR016
#SBATCH --time=07:00:00
#SBATCH --output=//mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC_trial/jobs_scripts/MIGEC-BulkTCR016.out
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

Rscript /home/moarias/scripts/github_submit/MIGEC/R/run_migec.R --ReportsDir=/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC_trial --ProjectID=BulkTCR016 --bcl2fastqSheet=/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/sequencing_data/04-17-2023/mkfastq/data/experiment_layout/NV103_sample_sheet_for_mkfastq.csv --bcl2fastqReports=/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/sequencing_data/04-17-2023/mkfastq/NV103/outs/fastq_path/BulkTCR016 --BarcodesSheet=/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC/metadata/BulkTCR016_barcodes_sheet_migec.csv --SampleSheet=/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/COVID-19/paper_developments/COVID-Vaccine-Assesment-Upper-Track/MIGEC/metadata/BulkTCR016_sample_sheet_migec.csv --TimeLimit='06:00:00' --Species=HomoSapiens --FileType=paired --Mask='1:1' --Quality='25,30' --MIGsize=2
