source("~/Documents/r_scripts_and_functions/r_functions/getting_cohort_info.R")
source("~/Documents/r_scripts_and_functions/r_functions/getting_pca_stats.R")

con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
con <- opencgaLogin(opencga = con, userid = "jambrose", passwd = getPass())

ldp_codes <- read_csv("~/Documents/trial_counting/trial.counting/mapping_plots/gmc_ldp_locations_and_codes.csv")

study <- 1000000038

cohort_details <- getting_cohort_info(con, study, upload_report)
pca_stats <- getting_pca_stats(con, study)

list <- read_csv("~/Documents/checking_bertha_stats/qms_jobs_for_john_2018-06-25.csv", col_names = c("tumour_sample_id", "qms_id"))

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

write_tsv(pca_stats_with_participant_id, paste0("~/Documents/checking_bertha_stats/pca_stats_",Sys.Date(),".tsv"))
