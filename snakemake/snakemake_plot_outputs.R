library(dplyr)
library(readr)
library(ggplot2)

# import arcadia font ---------------------------------------------------------

# this section of code defines an "Arcadia" theme layer for ggplot2 plots.
# It's a combination of theme_classic() and a definition for a brand-compatible font.
# If the font file is not available, theme_arcadia is synonymous with theme_classic().
# Additional theme definitions can be layered onto plots as needed.
# The font file is not defined as a snakemake input so that the snakefile will still execute in its absence.

if(file.exists("inputs/SuisseIntl-Regular.otf")){
  library(showtext)
  font_add("SuisseIntl-Regular", "inputs/SuisseIntl-Regular.otf")
  showtext::showtext_auto()
  theme_arcadia <- theme_classic() +
    theme(text = element_text(family = "SuisseIntl-Regular"))
} else {
  theme_arcadia <- theme_classic()
}

# read in and format data -------------------------------------------------

all <- read_tsv(snakemake@input[['all_outputs']])

# filter to genbank protein accessions that had uniprot accessions
all_with_class <- all %>%
  filter(!is.na(class)) %>%
  # edit class names so they'll be pretty on the plot
  mutate(class = gsub("arp", "actin related", class)) %>%
  mutate(class = gsub("actin like", "actin-like", class),
         class = gsub("actin related", "actin-related", class)) %>%
  # filter only to classes of interest
  filter(class %in% c("actin", "actin-like", "actin-related")) %>%
  # arrange to layer points so actin-like and actin-related can be seen
  arrange(class)

# figure 2: global sequence conservation ---------------------------------------

fig2b <- ggplot(all, aes(x = reorder(protein, avg_pid), y = avg_pid)) +
  geom_col(fill = "#5088C5", width = 2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Actin (AU)", y = "Average global\nsequence identity (%)") +
  theme_arcadia +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

pdf(snakemake@output[['fig2b']], width = 5, height = 3)
fig2b
dev.off()

fig2c <- ggplot(all, aes(x = avg_pid)) +
  geom_histogram(fill = "#5088C5", bins = 70) +
  scale_y_continuous(limits = c(0, 2500), expand = c(0, 0)) +
  labs(y = "Number of values", x = "Average global sequence identity (%)") +
  theme_arcadia

pdf(snakemake@output[['fig2c']], width = 5, height = 3)
fig2c
dev.off()

# figure 3: actin structure  ---------------------------------------------------

fig3e <- ggplot(all %>% filter(!is.na(evalue_transform)), 
                aes(x = reorder(protein, evalue_transform), y = evalue_transform)) +
  geom_col(fill = "#F28360", width=2) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  labs(x = "Actin (AU)", y = "Structural similarity") +
  theme_arcadia +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

pdf(snakemake@output[['fig3e']], width = 5, height = 3)
fig3e
dev.off()

fig3f <- ggplot(all %>% filter(!is.na(evalue_transform)), 
                aes(x = avg_pid, y = evalue_transform)) +
  geom_point(size = .5) +
  xlim(0, 100) + 
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  labs(x = "Average global percent identity (%)", y = "Structural similarity") +
  theme_arcadia

pdf(snakemake@output[['fig3f']], width = 5, height = 3)
fig3f
dev.off()

# figure 4: actin contacts -----------------------------------------------------

fig4c <- ggplot(all %>%
                  filter(!is.na(lat_fraction_matching)),
                aes(x = reorder(protein, lat_fraction_matching), y = lat_fraction_matching * 100)) +
  geom_col(fill = "#3B9886", width = 2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Actin (AU)", y = "Lateral contact\nconservation (%)") +
  theme_arcadia +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

pdf(snakemake@output[['fig4c']], width = 5, height = 3)
fig4c
dev.off()

fig4d <- ggplot(all %>%
                  filter(!is.na(lon_fraction_matching)), 
                aes(x = reorder(protein, lon_fraction_matching), y = lon_fraction_matching * 100)) +
  geom_col(fill = "#3B9886", width = 2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Actin (AU)", y = "Longitudinal contact\nconservation (%)") +
  theme_arcadia +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

pdf(snakemake@output[['fig4d']], width = 5, height = 3)
fig4d
dev.off()

fig4e <- ggplot(all %>%
                  filter(!is.na(w_avg_contacts)), 
                aes(x = reorder(protein, w_avg_contacts), y = w_avg_contacts * 100)) +
  geom_col(fill = "#3B9886", width = 2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Actin (AU)", y = "Total polymerization contact\nconservation (%)") +
  theme_arcadia +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

pdf(snakemake@output[['fig4e']], width = 5, height = 3)
fig4e
dev.off()

fig4f <- ggplot(all %>%
                  filter(!is.na(w_avg_contacts)), 
                aes(x = w_avg_contacts * 100)) +
  geom_histogram(fill = "#3B9886", bins = 18, position = position_dodge(width= .1 )) +
  scale_y_continuous(limits = c(0, 5000), expand = c(0, 0)) +
  xlim(10, 105) +
  labs(y = "Number of values", x = "Total polymerization contact conservation (%)") +
  theme_arcadia

pdf(snakemake@output[['fig4f']], width = 5, height = 3)
fig4f
dev.off()

fig4g <- ggplot(all %>% filter(!is.na(w_avg_contacts)),
                aes(x = avg_pid, y = w_avg_contacts * 100)) +
  geom_point(size = 1) +
  labs(x = "Average global sequence identity (%)", y = "Total polymerization contact\nconservation (%)") +
  theme_arcadia +
  scale_y_continuous(limits = c(0, 101), expand = c(0, 0))

pdf(snakemake@output[['fig4g']], width = 5, height = 3)
fig4g
dev.off()

fig4h <- ggplot(all %>% filter(!is.na(w_avg_contacts)),
                aes(x = evalue_transform, y = w_avg_contacts * 100)) +
  geom_point(size = 1) +
  labs(x = "Structural similarity", y = "Total polymerization contact\nconservation (%)") +
  theme_arcadia +
  scale_y_continuous(limits = c(0, 101), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 80))

pdf(snakemake@output[['fig4h']], width = 5, height = 3)
fig4h
dev.off()

# figure 5: actin atp ----------------------------------------------------------

fig5c <- ggplot(all %>%
                  filter(!is.na(atp_fraction_matching)),
                aes(x = reorder(protein, atp_fraction_matching), y = atp_fraction_matching * 100)) +
  geom_col(fill = "#3B9886", width = 2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Actin (AU)", y = "ATP-binding site\nconservation (%)") +
  theme_arcadia +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

pdf(snakemake@output[['fig5c']], width = 5, height = 3)
fig5c
dev.off()

fig5d <- ggplot(all %>%
                  filter(!is.na(atp_fraction_matching)), 
                aes(x = atp_fraction_matching * 100)) +
  geom_histogram(fill = "#3B9886", binwidth = 10, position = position_dodge(width= .1 )) +
  scale_y_continuous(limits = c(0, 10000), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(y = "Number of values", x = "ATP-binding site conservation (%)") +
  theme_arcadia

pdf(snakemake@output[['fig5d']], width = 5, height = 3)
fig5d
dev.off()

fig5e <- ggplot(all %>% filter(!is.na(atp_fraction_matching)),
                aes(x = avg_pid, y = atp_fraction_matching * 100)) +
  geom_point(size = 1) +
  labs(x = "Average global sequence identity (%)", y = "ATP-binding site\nconservation (%)") +
  theme_arcadia +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 105))

pdf(snakemake@output[['fig5e']], width = 5, height = 3)
fig5e
dev.off()

fig5f <- ggplot(all %>% filter(!is.na(atp_fraction_matching)),
                aes(x = evalue_transform, y = atp_fraction_matching * 100)) +
  geom_point(size = 1) +
  labs(x = "Structural similarity", y = "ATP-binding site\nconservation (%)") +
  theme_arcadia +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 80))

pdf(snakemake@output[['fig5f']], width = 5, height = 3)
fig5f
dev.off()

# figure 6: colored by actin type ----------------------------------------------

fig6a <- ggplot(all_with_class, aes(x = avg_pid, fill = class)) +
  geom_histogram(bins = 70, position = position_stack(reverse = TRUE)) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0, 0)) +
  labs(y = "Number of values", x = "Average global sequence identity (%)") +
  theme_arcadia + 
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.85)) +
  scale_fill_manual(values = c("#292928", "#9988DA", "#F7B846"),
                    labels = c('"Actin"', '"Actin-like"', '"Actin-related"'))

pdf(snakemake@output[['fig6a']], width = 5, height = 3)
fig6a
dev.off()

fig6b <-  ggplot(all_with_class, aes(x = avg_pid, y = evalue_transform, color = class)) +
  geom_point(size = .5) +
  xlim(0, 100) + 
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  labs(x = "Average global percent identity (%)", y = "Structural similarity") +
  theme_arcadia +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#292928", "#9988DA", "#F7B846"),
                     labels = c('"Actin"', '"Actin-like"', '"Actin-related"'))
# keep this to record how to add marginal density curves on the x and y axis on the scatter plot
# these are nice to illustrate where points are when there are two many points to discern each position
# they aren't in the pub at the moment because the arcadia font wouldn't work with gridExtra (a dependency of ggExtra)
#tmp <- ggExtra::ggMarginal(fig6b, type = "density", groupFill = T, groupColour = T)

pdf(snakemake@output[['fig6b']], width = 5, height = 3)
fig6b
dev.off()

fig6c <- ggplot(all_with_class %>% filter(!is.na(w_avg_contacts)),
                aes(x = avg_pid, y = w_avg_contacts * 100, color = class)) +
  geom_point(size = 1) +
  labs(x = "Average global sequence identity (%)", y = "Total polymerization contact\nconservation (%)") +
  theme_arcadia +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#292928", "#9988DA", "#F7B846"),
                     labels = c('"Actin"', '"Actin-like"', '"Actin-related"')) +
  scale_y_continuous(limits = c(0, 101), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(limits = c(0, 101), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) 

pdf(snakemake@output[['fig6c']], width = 5, height = 3)
fig6c
dev.off()

fig6d <- ggplot(all_with_class %>% filter(!is.na(atp_fraction_matching)),
                aes(x = avg_pid, y = atp_fraction_matching * 100, color = class)) +
  geom_point(size = 1) +
  labs(x = "Average global sequence identity (%)", y = "ATP-binding site\nconservation (%)") +
  theme_arcadia +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#292928", "#9988DA", "#F7B846"),
                     labels = c('"Actin"', '"Actin-like"', '"Actin-related"')) +
  scale_y_continuous(limits = c(0, 101), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(limits = c(0, 101), breaks = c(0, 25, 50, 75, 100))

pdf(snakemake@output[['fig6d']], width = 5, height = 3)
fig6d
dev.off()
