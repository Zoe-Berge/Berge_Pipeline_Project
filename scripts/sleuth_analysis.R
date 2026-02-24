library(sleuth)
library(dplyr)

sample_dirs <- list.dirs("results", recursive = FALSE, full.names = TRUE) #looks for the srr numbers in my results folder 
sample_dirs <- sample_dirs[file.exists(file.path(sample_dirs, "abundance.tsv"))] #keeps folders that have abundance.tsv file  
sample_names <- basename(sample_dirs)

if (length(sample_names) > 1) { #if there are samples, keeps my code from crashing  

    conditions <- ifelse(grepl("33|45", sample_names), "6dpi", "2dpi") #looks at the numbers and if they have a 33 or 45 they are 6dpi

    s2c <- data.frame( #dataframe telling sleuth where/what everything is
        sample = sample_names,
        condition = conditions,
        path = sample_dirs,
        stringsAsFactors = FALSE
    )

    #Sleuth Processing, liklihood ratio test  
    so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
    so <- sleuth_fit(so, ~condition, 'full')
    so <- sleuth_fit(so, ~1, 'reduced')
    so <- sleuth_lrt(so, 'reduced', 'full')

    sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

    #sleuth perfoming all of the tests with a q value of below 0.05
    sleuth_significant <- sleuth_table %>%
        dplyr::filter(qval < 0.05) %>%
        dplyr::select(target_id, test_stat, pval, qval) %>%
        dplyr::arrange(pval)
    
    #puts information in a table 
    write.table(sleuth_significant, file="results/significant_transcripts.temp", 
                sep="\t", quote=FALSE, row.names=FALSE)

} else {
    #if there is only one sample 
    write("target_id\ttest_stat\tpval\tqval\nTest_Mode\t0\t0\t0", 
          file="results/significant_transcripts.temp")
}


