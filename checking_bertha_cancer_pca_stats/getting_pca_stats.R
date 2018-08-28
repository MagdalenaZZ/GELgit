getting_pca_stats <- function(con, study){
  require(janitor)
  require(getPass)
  require(tidyverse)
  require(opencgaR)

  file_stats_tumour <- fileClient(OpencgaR = con, action = "search", params = list(study = study, format = "BAM", include = "name,stats.SAMPLE_METADATA.sample_library_type,stats.AT_GC_DROP.gc_drop,stats.AT_GC_DROP.at_drop,stats.WHOLE_GENOME_COVERAGE.coverageSummary,stats.SAMTOOLS_STATS_CALCULATIONS_ALL.SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES_PERCENT,stats.SAMTOOLS_STATS_ALL.SAMTOOLS_INSERT_SIZE_AVERAGE,stats.SAMTOOLS_STATS_CALCULATIONS_ALL.SAMTOOLS_READS_MAPPED_PERCENT,stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SNVS,stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_INDELS,stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SNVS", status = "READY")) %>% 
    select(id, name, stats) %>%
    jsonlite::flatten()
  
  pca_stats <- tibble(tumour_bam_id = file_stats_tumour$name,
                      library_type = file_stats_tumour$stats.SAMPLE_METADATA.sample_library_type,
                      gc_drop = file_stats_tumour$stats.AT_GC_DROP.gc_drop,
                      at_drop = file_stats_tumour$stats.AT_GC_DROP.at_drop,
                      coverage_homogeneity = map_dbl(file_stats_tumour$stats.WHOLE_GENOME_COVERAGE.coverageSummary, ~ ifelse(!is.null(.x), .x[which(.x$scope == "autosomes"),"localRMSD"], NA)),
                      chimeric_percentage = file_stats_tumour$stats.SAMTOOLS_STATS_CALCULATIONS_ALL.SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES_PERCENT,
                      average_fragment_size = file_stats_tumour$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_INSERT_SIZE_AVERAGE,
                      mapping_rate = file_stats_tumour$stats.SAMTOOLS_STATS_CALCULATIONS_ALL.SAMTOOLS_READS_MAPPED_PERCENT,
                      somatic_snvs = file_stats_tumour$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SNVS,
                      somatic_indels = file_stats_tumour$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_INDELS) %>%
    mutate(tumour_sample_id = str_replace_all(tumour_bam_id, ".bam", ""),
           somatic_snvs_log10 = round(log(somatic_snvs, 10), digits = 2),
           somatic_indels_log10 = round(log(somatic_indels, 10), digits = 2),
           snvs_per_mb_log10 = round(log(somatic_snvs/(3077073773/1000000), 10), digits = 2)) %>% 
    select(-c(tumour_bam_id, somatic_snvs, somatic_indels)) %>% 
    select(tumour_sample_id, everything())
  
  return(pca_stats)
}