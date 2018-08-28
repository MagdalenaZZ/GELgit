getting_cohort_info <- function(con, study, upload_report){
  require(janitor)
  require(getPass)
  require(opencgaR)
  require(tidyverse)
  
  cohorts <- cohortClient(OpencgaR = con, action = "search", params = list(study = study))
  samples <- sampleClient(con, action = "search", params = list(study = study))
  
  samples_names <- samples %>% 
    select(id, name) %>% 
    rename(sample_id = id, sample_name = name)
  
  samples_attributes <- samples %>% 
    select(id, attributes) %>%
    jsonlite::flatten() %>% 
    select(id, attributes.processes)
  
  #gets the latest cohort by the final part of the cohort name - as advised by Kenan
  cohorts_and_info <- cohorts %>% 
    select(id, name, annotationSets) %>% 
    unnest() %>% 
    select(id, name, annotations) %>%
    rename(cohort_id = id, cohort_name = name) %>% 
    unnest() %>% 
    filter(name == "individualId" | 
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

  cohorts_and_samples <- cohorts %>% 
    select(id, name, samples) %>% 
    unnest() %>% 
    select(-c(release, version, somatic)) %>% 
    unique() %>% 
    rename(cohort_id = id, cohort_name = name, sample_id = id1) %>% 
    as.tibble()
  
  cancer_delivery_ids <- upload_report %>% 
    select(Platekey, DeliveryID, Type, `Delivery Date`, `Delivery Version`) %>% 
    rename(sample_name = Platekey, delivery_id = DeliveryID, sample_type = Type, delivery_date = `Delivery Date`, delivery_version = `Delivery Version`) %>% 
    filter(sample_type != "rare disease") %>% 
    group_by(sample_name) %>% 
    filter(delivery_date == max(delivery_date),
           delivery_version == max(delivery_version)) %>%
    ungroup() %>% 
    mutate(sample_type = str_replace_all(sample_type, "cancer ", ""))
  
  bertha_versions <- map2_df(samples_attributes$id, samples_attributes$attributes.processes,
                             ~ tibble(id = .x,
                                      process_name = ifelse(!is.null(.y), .y$name[order(.y$processDate, decreasing = TRUE)], NA),
                                      process_version = ifelse(!is.null(.y), .y$version[order(.y$processDate, decreasing = TRUE)], NA),
                                      process_date = ifelse(!is.null(.y), .y$processDate[order(.y$processDate, decreasing = TRUE)], NA)))
  
  cohort_info <- left_join(cohorts_and_info, cohorts_and_samples, by = c("cohort_id" = "cohort_id")) %>% 
    left_join(., samples_names, by = c("sample_id" = "sample_id")) %>% 
    left_join(., cancer_delivery_ids, by = c("sample_name" = "sample_name")) %>% 
    left_join(., bertha_versions, by = c("sample_id" = "id")) %>% 
    filter(readyForAnalysis == TRUE) %>% 
    select(-readyForAnalysis) %>% 
    rename(ldp = center,
           individual_id = individualId,
           primary_diagnosis_disease = primaryDiagnosisDisease,
           primary_diagnosis_subdisease = primaryDiagnosisSubDisease,
           year_of_birth = yearOfBirth)
  
  return(cohort_info)
}

