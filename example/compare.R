library(ggplot2)
library(tidyr)
library(dplyr)
library(skimr)

cont_1000g_dv <- read.csv("1000G_DeepVar.sce.tsv", header=T, sep="\t")
cont_1000g_GATK <- read.csv("1000G_GATK.sce.tsv", header=T, sep="\t")

cont_1000g_dv$cohort <- "1000G_DV"
cont_1000g_GATK$cohort <- "1000G_GATK"

df <- rbind(cont_1000g_dv, cont_1000g_GATK)
df_plot <- df %>% rename( SAMPLE = X.SAMPLE ) %>% select(-HQ_HOM_RATE,-HQ_HET_RATE) %>%
  pivot_longer(HQ_HOM:INCONSISTENT_AB_HET_RATE, names_to="metric", values_to = "value")

p <- ggplot(df_plot, aes(x=value, color=cohort)) + geom_density() + facet_wrap(~metric, scales="free") + theme_bw()
ggsave(plot = p, filename = "Metric_distribution.pdf", device="pdf", height = 7, width = 10)

sum_df <- df %>% group_by(cohort) %>% skim(HQ_HOM:INCONSISTENT_AB_HET_RATE)
write.table(sum_df, file = "Summary_metric_distribution.tsv", sep="\t", row.names=F, quote=F)
