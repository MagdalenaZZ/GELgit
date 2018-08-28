
require(opencgaR)
require(getPass)

require(janitor)
require(tidyverse)


cases <- c("215999003","215999008","220000517","220000521","220000522","220000576","220000615","220000647","220000663","220000730","220000743","220000811","220000884","220000885","220000955","220000962","220000980","220000986","220000996")
casescases <-  c("LP3000396-DNA_A01","LP3000396-DNA_A02","LP3000396-DNA_B01","LP3000396-DNA_B02","LP3000396-DNA_C01","LP3000396-DNA_C02","LP3000396-DNA_D01","LP3000396-DNA_D02","LP3000396-DNA_E01","LP3000396-DNA_E02","LP3000396-DNA_F01","LP3000396-DNA_F02","LP3000396-DNA_G01","LP3000396-DNA_G02","LP3000396-DNA_H01","LP3000396-DNA_H02")


con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
con <- opencgaLogin(opencga = con, userid = "mzarowiecki", passwd = getPass())
study <- 1000000038

#upload_report <- "/Users/magz/git/checking_bertha_cancer_pca_stats/upload_report_latest.txt"
#source("~/Documents/r_scripts_and_functions/r_functions/getting_cohort_info.R")
#source("~/Documents/r_scripts_and_functions/r_functions/getting_pca_stats.R")


cohorts <- cohortClient(OpencgaR = con, action = "search", params = list(study = study))
#samples <- sampleClient(con, action = "search", params = list(study = study))


#gets the latest cohort by the final part of the cohort name - as advised by Kenan
cohorts_and_info2 <- cohorts %>% 
	select(id, name, annotationSets) %>% 
    	unnest() %>% 
    	select(id, name, annotations) %>%
    	rename(cohort_id = id, cohort_name = name) %>% 
    	unnest() %>% 
    	filter(name == "individualId" | 
             name == "matchedSamples" | 
             name == "primaryDiagnosisDisease" | 
             name == "primaryDiagnosisSubDisease" | 
             name == "readyForAnalysis" | 
             name == "center" | 
             name == "sex" | 
             name == "yearOfBirth") %>% 
    	spread(name, value) %>% 
	as.tibble() %>% 
    	mutate_if(is.list, as.character) %>% 
    	mutate(readyForAnalysis = as.logical(readyForAnalysis),
		sex = as.factor(sex),
           	yearOfBirth = as.integer(yearOfBirth),
           	cohort_version = str_extract(cohort_name, "(?<=_)[^_]*$")) %>% 
    	group_by(individualId) %>% 
    	filter(cohort_version == max(cohort_version)) %>%
    	ungroup() %>% 
    	select(-c(cohort_name, cohort_version)) %>% 
	unique()



cohorts_and_info2$germlineSampleId <- sapply(strsplit(as.character(cohorts_and_info2$matchedSamples), "\""), `[`, 2)
cohorts_and_info2$tumourSampleId <- sapply(strsplit(as.character(cohorts_and_info2$matchedSamples), "\""), `[`, 4)



res <- list()
j=1


# Pick only samples in cases
cohorts_and_info2[cohorts_and_info2$tumourSampleId  %in% cases,c(10,11,5,6)]




