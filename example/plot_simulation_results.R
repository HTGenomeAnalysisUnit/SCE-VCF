library(ggplot2)
library(tidyr)
library(dplyr)

df <- read.table("simulation_results.tsv", header=T, comment.char = "", sep="\t")
df <- df %>% separate(X.SAMPLE, into=c("id","original","cont_fraction"),sep="_") %>%
  select(-original) %>% replace_na(list(cont_fraction="0"))
df <- df %>% mutate(across(cont_fraction:INCONSISTENT_AB_HET_RATE, as.numeric))

df_plot <- df %>% select(cont_fraction,CHARR, HETEROZYGOSITY_RATE, INCONSISTENT_AB_HET_RATE) %>% 
  pivot_longer(CHARR:INCONSISTENT_AB_HET_RATE, names_to = "metric", values_to = "value")
p <- ggplot(df_plot, aes(x=cont_fraction, y=value)) + geom_point() + geom_smooth() + 
    facet_wrap(~metric, scales="free_y") + theme_bw() 
ggsave(p, file="simulation_metric_plot.pdf", device="pdf", height=4, width=8, dpi=150)
ggsave(p, file="simulation_metric_plot.png", device="png", height=4, width=8, dpi=150)
