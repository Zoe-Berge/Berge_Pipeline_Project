library(sleuth)
library(dplyr)

sample_dirs <- list.dirs("results", recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[file.exists(file.path(sample_dirs, "abundance.tsv"))]
sample_names <- basename(sample_dirs)

if (length(sample_names) > 1) { #if there are samples 

    conditions <- ifelse(grepl("33|45", sample_names), "6dpi", "2dpi")

    s2c <- data.frame(
        sample = sample_names,
        condition = conditions,
        path = sample_dirs,
        stringsAsFactors = FALSE
    )

    #Sleuth Processing
    so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
    so <- sleuth_fit(so, ~condition, 'full')
    so <- sleuth_fit(so, ~1, 'reduced')
    so <- sleuth_lrt(so, 'reduced', 'full')

    sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

    sleuth_significant <- sleuth_table %>%
        dplyr::filter(qval < 0.05) %>%
        dplyr::select(target_id, test_stat, pval, qval) %>%
        dplyr::arrange(pval)

    write.table(sleuth_significant, file="results/significant_transcripts.temp", 
                sep="\t", quote=FALSE, row.names=FALSE)

} else {
    write("target_id\ttest_stat\tpval\tqval\nTest_Mode\t0\t0\t0", 
          file="results/significant_transcripts.temp")
}


