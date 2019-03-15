library(tidyverse)

freq = c(
  327178,-854269,
  26511,-514729,
  36945,-893219,
  100624,-593706,
  4753,-325979,
  56037,-650118,
  98178,-1016677,
  125652,-966034
) / 52084138

target_sites_df = tibble(
  tissue = c(
    'Cerebellum','Cerebellum',
    'Medulla Oblongata','Medulla Oblongata',
    'Frontal lobe','Frontal lobe',
    'Parietal lobe','Parietal lobe',
    'Heart','Heart',
    'Liver','Liver',
    'Kidney','Kidney',
    'HEK293','HEK293'
  ),
  gain_or_loss = c(
    'Gain','Loss',
    'Gain','Loss',
    'Gain','Loss',
    'Gain','Loss',
    'Gain','Loss',
    'Gain','Loss',
    'Gain','Loss',
    'Gain','Loss'
  ),
  freq = c(
    0.628172055,-1.640171140,
    0.050900334,-0.988264412,
    0.070933304,-1.714953985,
    0.193195095,-1.139897909,
    0.009125619,-0.625870011,
    0.107589378,-1.248207276,
    0.188498848,-1.951989683,
    0.241248113,-1.854756625
  )
)

target_sites_df$tissue = as.factor(target_sites_df$tissue)

p <- 
  target_sites_df %>% 
  ggplot(aes(x = tissue, y = freq, group = gain_or_loss, fill = gain_or_loss)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_x_discrete(limits = target_sites_df$tissue) +
  # another trick!
  scale_y_continuous(
    breaks = seq(-2,2,0.4),
    labels = abs(seq(-2,2,0.4)),
    limits = c(-2,2)
  ) +
  labs(x = "Sample", 
       y = "Percentage of total original sites", 
       title = "Total miRNA target site gain and loss",
       subtitle = "Sequencing Depth: 20 million reads"
  ) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill =  "grey90")) +
  # reverse the order of items in legend
  # guides(fill = guide_legend(reverse = TRUE)) +
  # change the default colors of bars
  scale_fill_manual(values=c(
    rgb(0,0,1,0.55),
    rgb(1,0,0,0.55)
  ),
  name="",
  breaks=c("Gain", "Loss"),
  labels=c("Gain", "Loss")) + theme_classic()

ggsave('results/plots/figure_4.png',plot=p, scale=0.75)
