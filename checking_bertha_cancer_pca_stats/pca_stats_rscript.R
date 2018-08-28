#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Four arguments must be supplied (i. latest upload_report, ii. list of ldp codes, iii. tumour sample and QMS ticket input file and, iv. openCGA userid)\n", call. = FALSE)
} else if (length(args) == 4) {
  upload_report_path <- args[1]
  ldp_codes_filepath <- args[2]
  input_filename <- args[3]
  userid <- args[4]
}

library(getPass)
library(tidyverse)
library(opencgaR)
#THE BELOW FILEPATHS MUST BE CHANGED TO MATCH THE LOCATION OF THESE FUNCTION SCRIPTS ON THE LOCAL COMPUTER OF THE USER
source("~/Documents/r_scripts_and_functions/r_scripts/checking_bertha_stats/getting_cohort_info.R")
source("~/Documents/r_scripts_and_functions/r_scripts/checking_bertha_stats/getting_pca_stats.R")

output_filename <- paste0(str_replace_all(input_filename, ".tsv", ""), ".output.")

con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
con <- opencgaLogin(opencga = con, userid = userid, passwd = getPass())

#scp jambrose@10.5.8.15:/genomes/resources/upload_reports/upload_report.2018-06-13.txt ~/
upload_report <- read_tsv(upload_report_path, skip = 15)

ldp_codes <- read_csv(ldp_codes_filepath)

study <- 1000000038

cohort_details <- getting_cohort_info(con, study, upload_report)
pca_stats <- getting_pca_stats(con, study)

list <- read_tsv(input_filename, col_names = c("tumour_sample_id", "qms_id"))

pca_stats_with_participant_id <- left_join(list, pca_stats, by = c("tumour_sample_id" = "tumour_sample_id")) %>% 
  left_join(., cohort_details, by = c("tumour_sample_id" = "sample_name")) %>%
  left_join(., ldp_codes, by = c("ldp" = "ODS Code")) %>% 
  mutate(library_type = str_replace_all(library_type, "TruSeq PCR-Free High Throughput", "PCR-Free"),
         library_type = str_replace_all(library_type, "pcrfree", "PCR-Free"),
         library_type = str_replace_all(library_type, "TruSeq Nano High Throughput", "Nano"),
         date = Sys.Date(),
         qms_ticket = paste0("https://jira.extge.co.uk/browse/",qms_id),
         qms_id = str_replace_all(qms_id, "QMS-", "")) %>% 
  rename(gmc_code = GMC_code, hospital = ldp) %>%
  select(tumour_sample_id,
         gmc_code,
         individual_id,
         sample_id,
         qms_ticket,
         qms_id,
         hospital, 
         everything(),
         -starts_with("cohort"),
         -starts_with("germline_sample"),
         -Site,
         -GMC,
         -starts_with("process"),
         -year_of_birth,
         -primary_diagnosis_subdisease,
         -sex,
         -starts_with("delivery"),
         -sample_type,
         -snvs_per_mb_log10) %>% 
  arrange(gmc_code)

head(pca_stats_with_participant_id)

write_tsv(pca_stats_with_participant_id, paste0(output_filename, Sys.Date(), ".tsv"))
