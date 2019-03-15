library(tidyverse)

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
    327178,-854269,
    26511,-514729,
    36945,-893219,
    100624,-593706,
    4753,-325979,
    56037,-650118,
    98178,-1016677,
    125652,-966034
    ) / 1000000
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
    breaks = seq(-1.2,1.2,0.3),
    labels = abs(seq(-1.2,1.2,0.3)),
    limits = c(-1.2,1.2)
  ) +
  labs(x = "Sample", 
       y = "Count (millions)", 
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

ggsave('results/plots/sf_1.png', plot = p, scale = 0.75)

