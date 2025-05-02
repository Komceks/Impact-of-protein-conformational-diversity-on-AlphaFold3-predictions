# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggpubr)
library(gridExtra)
# Load the RMSD data
data <- read.csv("RMSD_results.csv", sep=',')  # Replace with your actual file path

# Add classification for A and B
data <- data %>%
  mutate(
    Closer_To_A = ifelse(HOLO_vs_APO_ALPHA < APO_vs_APO_ALPHA, "Holo", "Apo"),
    Closer_To_B = ifelse(HOLO_vs_HOLO_ALPHA < APO_vs_HOLO_ALPHA, "Holo", "Apo")
  )

figure0 <- ggplot(data, aes(x = APO_vs_HOLO)) +
  geom_density(alpha = 0.3, color = "black", fill = "black", adjust = 1) +
  ggtitle("RCSB struktūrų Cα-RMSD pasiskirstymas") +
  xlab(expression(RMSD[Apo~vs~Holo]~(Å))) +
  ylab("Tankis") +
  scale_x_continuous(breaks = seq(0,70,2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.6),
                     breaks = seq(0, 0.6, 0.05),
                     expand = c(0,0)) +
  coord_cartesian(xlim = c(0,42)) +      # <-- zoom without dropping data
  theme_classic(base_size = 13) +
  theme(plot.margin = unit(c(.2,.2,.2,.2), "cm"))

print(figure0)
summary(data$APO_vs_HOLO)
print(sum(data$APO_vs_HOLO < 2) / nrow(data))
ggsave("Stats/fig0_density_a3_org.png", plot = figure0, width = 8, height = 6)

d <- density(data$APO_ALPHA_vs_HOLO_ALPHA, from = 0, to = 20, adjust = 1)

# 2. turn it into a data.frame
df_d <- data.frame(x = d$x, y = d$y)

# 3. plot with geom_area()
figure0_2 <- ggplot(df_d, aes(x = x, y = y)) +
  geom_area(fill = "black", color = "black", alpha = 0.3) +
  ggtitle("Prognozuotų struktūrų Cα-RMSD pasiskirstymas") +
  xlab(expression(RMSD[Apo_alpha~vs~Holo_alpha]~(Å))) +
  ylab("Tankis") +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, 2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 2),
    breaks = seq(0, 2, 0.1),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 13) +
  theme(plot.margin = unit(c(.2, .2, .2, .2), "cm"))

print(figure0_2)
summary(data$APO_ALPHA_vs_HOLO_ALPHA)
print(sum(data$APO_ALPHA_vs_HOLO_ALPHA < 2) / nrow(data))
ggsave("Stats/fig0_density_a3_alpha.png", plot = figure0_2, width = 8, height = 6)

combined_figure = grid.arrange(
  figure0, figure0_2,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig0_density_a3_combined.png", plot = combined_figure, width = 16, height = 8)


plddt_data <- read.csv("plddt_results.tsv", sep = "\t") %>%
  mutate(plddt_mean = rowMeans(across(PLDDT_0:PLDDT_4), na.rm = TRUE))

# 2. Now summarise over that new column
plddt_summary <- plddt_data %>%
  summarise(
    overall_mean     = mean(plddt_mean, na.rm = TRUE),
    total_rows       = n(),
    count_above85    = sum(plddt_mean >= 85, na.rm = TRUE),
    percent_above85  = 100 * mean(plddt_mean >= 85, na.rm = TRUE)
  )

# 3. Inspect
print(plddt_summary)
# FIG 1 - APO_ALPHA TO APO AND HOLO

figure1a_full <- data %>% 
  ggplot(aes(x = APO_vs_APO_ALPHA, y = HOLO_vs_APO_ALPHA, color = Closer_To_A)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = c("Holo" = "blue", "Apo" = "red")) +
  labs(x = expression(RMSD[~vs~Apo]~(Å)),
       y = expression(RMSD[~vs~Holo]~(Å))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
  scale_x_continuous(breaks = seq(0, 50, 2), limits = c(0, 50), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 50, 2), limits = c(0, 50), expand = c(0, 0)) +
  theme_minimal() +
  ggtitle('1A: prognozuotų apo struktūrų panašumo įvertinimas su apo ir holo \n eksperimentinėmis struktūromis (pilnas)')

figure1b_zoom <- data %>% 
  ggplot(aes(x = APO_vs_APO_ALPHA, y = HOLO_vs_APO_ALPHA, color = Closer_To_A)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Holo" = "blue", "Apo" = "red")) +
  labs(x = expression(RMSD[~vs~Apo]~(Å)),
       y = expression(RMSD[~vs~Holo]~(Å))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
  coord_cartesian(xlim = c(0,3.5), ylim = c(0,3.5)) +
  scale_x_continuous(breaks = seq(0, 3.5, 0.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 3.5, 0.5), expand = c(0, 0)) +
  theme_minimal() +
  ggtitle('1B: prognozuotų apo struktūrų panašumo įvertinimas su apo ir holo \n eksperimentinėmis struktūromis (priartinta iki 0–3.5 Å)')

print(figure1a_full)
print(figure1b_zoom)

ggsave("Stats/fig1a_a3_full.png", plot = figure1a_full, width = 16, height = 8)
ggsave("Stats/fig1b_a3_zoom.png", plot = figure1b_zoom, width = 16, height = 8)

combined_figure = grid.arrange(
  figure1a_full, figure1b_zoom,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig1_a3_combined.png", plot = combined_figure, width = 16, height = 8)

table(data$Closer_To_A)

# # FIG 2 - HOLO_ALPHA TO APO AND HOLO
figure2a_full <- data %>% 
  ggplot(aes(x = APO_vs_HOLO_ALPHA, y = HOLO_vs_HOLO_ALPHA, color = Closer_To_B)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = c("Holo" = "blue", "Apo" = "red")) +
  labs(x = expression(RMSD[~vs~Apo]~(Å)),
       y = expression(RMSD[~vs~Holo]~(Å))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
  scale_x_continuous(breaks = seq(0, 50, 2), limits = c(0, 50), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 50, 2), limits = c(0, 50), expand = c(0, 0)) +
  theme_minimal() +
  ggtitle('2A: prognozuotų holo struktūrų panašumo įvertinimas su apo ir holo \n eksperimentinėmis struktūromis (pilnas)')

figure2b_zoom <- data %>% 
  ggplot(aes(x = APO_vs_HOLO_ALPHA, y = HOLO_vs_HOLO_ALPHA, color = Closer_To_B)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Holo" = "blue", "Apo" = "red")) +
  labs(x = expression(RMSD[~vs~Apo]~(Å)),
       y = expression(RMSD[~vs~Holo]~(Å))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
  coord_cartesian(xlim = c(0,3.5), ylim = c(0,3.5)) +
  scale_x_continuous(breaks = seq(0, 3.5, 0.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 3.5, 0.5), expand = c(0, 0)) +
  theme_minimal() +
  ggtitle('2B: prognozuotų holo struktūrų panašumo įvertinimas su apo ir holo \n eksperimentinėmis struktūromis (priartinta iki 0–3.5 Å)')

print(figure2a_full)
print(figure2b_zoom)

ggsave("Stats/fig2a_a3_full.png", plot = figure2a_full, width = 16, height = 8)
ggsave("Stats/fig2b_a3_zoom.png", plot = figure2b_zoom, width = 16, height = 8)

combined_figure = grid.arrange(
  figure2a_full, figure2b_zoom,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig2_a3_combined.png", plot = combined_figure, width = 16, height = 8)

table(data$Closer_To_B)

table(data$Closer_To_A == "Apo" & data$Closer_To_B == "Holo")
print(mean(data$HOLO_vs_APO_ALPHA))
print(mean(data$APO_vs_APO_ALPHA))
wilcoxon_result <- wilcox.test(
  c(data$APO_vs_APO_ALPHA), 
  c(data$HOLO_vs_APO_ALPHA),
  alternative = "less"
)
# W = 4809, p-value = 0.9851
cat("Wilcoxon Signed-Rank Test:\n")
print(wilcoxon_result)

print(mean(data$HOLO_vs_HOLO_ALPHA))
print(mean(data$APO_vs_HOLO_ALPHA))
wilcoxon_result <- wilcox.test(
  c(data$HOLO_vs_HOLO_ALPHA), 
  c(data$APO_vs_HOLO_ALPHA),
  alternative = "less"
)
# W = 2845.5, p-value = 0.000286
cat("Wilcoxon Signed-Rank Test:\n")
print(wilcoxon_result)

# FIG 4
# Reshape data to have PDB_ID and corresponding APO_vs_HOLO
reshaped_rmsd <- data %>%
  select(Apo_PDB_ID, Holo_PDB_ID, APO_vs_HOLO) %>% 
  mutate(APO_vs_HOLO = as.numeric(gsub(",", ".", APO_vs_HOLO))) %>%  # Ensure numeric format
  pivot_longer(
    cols = c(Apo_PDB_ID, Holo_PDB_ID),
    names_to = "PDB_Type",
    values_to = "PDB_ID"
  )

# Merge with plDDT data
merged_data <- reshaped_rmsd %>%
  # build a new PDB_NAME that is exactly "XXXX_APO" or "XXXX_HOLO"
  mutate(
    PDB_NAME = paste0(
      toupper(PDB_ID),
      "_",
      toupper(strtrim(PDB_Type, 4))
    )
  ) %>%
  inner_join(plddt_data, by = "PDB_NAME")

# Create the scatter plot
figure4a_full <- ggplot(merged_data, aes(x = APO_vs_HOLO, y = plddt_mean)) +
  geom_point(fill = 'black', color = 'black', alpha = 0.3, size = 4) +
  geom_smooth(method = "lm", color = "blue") +
  
  # Axis labels and titles
  labs(
    x = expression(RMSD[Apo~vs~Holo]~(Å)),
    y = expression(plDDT)
  ) +
  
  # Axis scales
  scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), expand = c(0, 0)) +
  
  # Theme customization
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = 'plain', size = 12, margin = margin(t = 10, r = 8, b = 0, l = 0)),
    axis.title.y = element_text(face = 'plain', size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    plot.margin = unit(c(.2, .2, .2, .2), "cm"),
    plot.title = element_text(face = 'plain', size = 15)
  ) +
  
  ggtitle("plDDT įverčio pokytis didėjant konformacijos įvairovei (pilnas)")

filtered_merged_data <- subset(merged_data, APO_vs_HOLO < 10)


figure4b_filtered <- ggplot(filtered_merged_data, aes(x = APO_vs_HOLO, y = plddt_mean)) +
  geom_point(fill = 'black', color = 'black', alpha = 0.3, size = 4) +
  geom_smooth(method = "lm", color = "blue") +
  
  # Axis labels and titles
  labs(
    x = expression(RMSD[Apo~vs~Holo]~(Å)),
    y = expression(plDDT)
  ) +
  
  # Axis scales
  scale_x_continuous(breaks = seq(0, 5, 0.5), limits = c(0, 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), expand = c(0, 0)) +
  
  # Theme customization
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = 'plain', size = 12, margin = margin(t = 10, r = 8, b = 0, l = 0)),
    axis.title.y = element_text(face = 'plain', size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    plot.margin = unit(c(.2, .2, .2, .2), "cm"),
    plot.title = element_text(face = 'plain', size = 15)
  ) +
  
  ggtitle("plDDT įverčio pokytis didėjant konformacijos įvairovei (iki 5 Å)")

print(figure4a_full)
ggsave("Stats/fig4a_full_a3_plDDT_RMSD.png", plot = figure4a_full, width = 8, height = 6)
print(figure4b_filtered)
ggsave("Stats/fig4b_filtered_a3_plDDT_RMSD.png", plot = figure4b_filtered, width = 8, height = 6)

combined_figure = grid.arrange(
  figure4a_full, figure4b_filtered,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig4a_combined.png", plot = combined_figure, width = 16, height = 8)


# p-value = 0.4145, cor = -0.08706841
cor.test(merged_data$APO_vs_HOLO, merged_data$plddt_mean,
         method = 'pearson')

filtered <- subset(merged_data, APO_vs_HOLO < 5)

# 2. Run Pearson correlation on the filtered subset

cor.test(
  filtered$APO_vs_HOLO,
  filtered$plddt_mean,
  method = "pearson"
)

data$Min_RMSD <- pmin(data$APO_vs_APO_ALPHA, data$HOLO_vs_HOLO_ALPHA, na.rm = TRUE)
data$Max_RMSD <- pmax(data$APO_vs_APO_ALPHA, data$HOLO_vs_HOLO_ALPHA, na.rm = TRUE)


# Create the scatter plot
figure5a_full_max <- ggplot(data, aes(x = APO_vs_HOLO, y = Max_RMSD)) +
  geom_point(fill = 'black', color = 'black', alpha = 0.3, size = 4) +
  geom_smooth(method = "lm", color = "blue") +
  
  # Axis labels and titles
  labs(
    x = expression(RMSD[Apo~vs~Holo]~(Å)),
    y = expression(Didžiausias~RMSD~įvertis)
  ) +
  
  # Axis scales
  scale_x_continuous(breaks = seq(0, 60, 2), limits = c(0, 60), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 20, 1), limits = c(0, 20), expand = c(0, 0)) +
  
  # Theme customization
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = 'plain', size = 12, margin = margin(t = 10, r = 8, b = 0, l = 0)),
    axis.title.y = element_text(face = 'plain', size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    plot.margin = unit(c(.2, .2, .2, .2), "cm"),
    plot.title = element_text(face = 'plain', size = 15)
  ) +
  
  ggtitle("A: didžiausio RMSD įverčio pokytis tarp prognozuotų ir eksperimentinių\n struktūrų didėjant struktūrų skirtumui tarp eksperimentinių baltymų struktūrų")


print(figure5a_full_max)
ggsave("Stats/fig5a_full_max_a3_RMSD_RMSD.png", plot = figure5a_full_max, width = 8, height = 6)


figure5b_full_min <- ggplot(data, aes(x = APO_vs_HOLO, y = Min_RMSD)) +
  geom_point(fill = 'black', color = 'black', alpha = 0.3, size = 4) +
  geom_smooth(method = "lm", color = "blue") +
  
  # Axis labels and titles
  labs(
    x = expression(RMSD[Apo~vs~Holo]~(Å)),
    y = expression(Mažiausias~RMSD~įvertis)
  ) +
  
  # Axis scales
  scale_x_continuous(breaks = seq(0, 60, 2), limits = c(0, 60), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 20, 1), limits = c(0, 20), expand = c(0, 0)) +
  
  # Theme customization
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = 'plain', size = 12, margin = margin(t = 10, r = 8, b = 0, l = 0)),
    axis.title.y = element_text(face = 'plain', size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    plot.margin = unit(c(.2, .2, .2, .2), "cm"),
    plot.title = element_text(face = 'plain', size = 15)
  ) +
  
  ggtitle("B: mažiausio RMSD įverčio pokytis tarp prognozuotų ir eksperimentinių\n struktūrų didėjant struktūrų skirtumui tarp eksperimentinių baltymų struktūrų")


print(figure5b_full_min)
ggsave("Stats/fig5b_full_min_a3_RMSD_RMSD.png", plot = figure5b_full_min, width = 8, height = 6)


# p-value = 0.0937, cor = 0.1777597
cor.test(data$APO_vs_HOLO, data$Max_RMSD,
         method = 'pearson')

# p-value = 0.9968, cor = 0.0004326738
cor.test(data$APO_vs_HOLO, data$Min_RMSD,
         method = 'pearson')

filtered_data <- subset(data, APO_vs_HOLO < 5)

# Create the scatter plot
figure5c_filtered_max <- ggplot(filtered_data, aes(x = APO_vs_HOLO, y = Max_RMSD)) +
  geom_point(fill = 'black', color = 'black', alpha = 0.3, size = 4) +
  geom_smooth(method = "lm", color = "blue") +
  
  # Axis labels and titles
  labs(
    x = expression(RMSD[Apo~vs~Holo]~(Å)),
    y = expression(Didžiausias~RMSD~įvertis)
  ) +
  
  # Axis scales
  scale_x_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0, 7), expand = c(0, 0)) +
  
  # Theme customization
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = 'plain', size = 12, margin = margin(t = 10, r = 8, b = 0, l = 0)),
    axis.title.y = element_text(face = 'plain', size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    plot.margin = unit(c(.2, .2, .2, .2), "cm"),
    plot.title = element_text(face = 'plain', size = 15)
  ) +
  
  ggtitle("C: didžiausio RMSD įverčio pokytis tarp prognozuotų ir eksperimentinių\n struktūrų didėjant struktūrų skirtumui tarp eksperimentinių baltymų struktūrų")


print(figure5c_filtered_max)
ggsave("Stats/fig5c_full_max_a3_RMSD_RMSD.png", plot = figure5c_filtered_max, width = 8, height = 6)

figure5d_filtered_min <- ggplot(filtered_data, aes(x = APO_vs_HOLO, y = Min_RMSD)) +
  geom_point(fill = 'black', color = 'black', alpha = 0.3, size = 4) +
  geom_smooth(method = "lm", color = "blue") +
  
  # Axis labels and titles
  labs(
    x = expression(RMSD[Apo~vs~Holo]~(Å)),
    y = expression(Mažiausias~RMSD~įvertis)
  ) +
  
  # Axis scales
  scale_x_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0, 7), expand = c(0, 0)) +
  
  # Theme customization
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = 'plain', size = 12, margin = margin(t = 10, r = 8, b = 0, l = 0)),
    axis.title.y = element_text(face = 'plain', size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    plot.margin = unit(c(.2, .2, .2, .2), "cm"),
    plot.title = element_text(face = 'plain', size = 15)
  ) +
  
  ggtitle("D: mažiausio RMSD įverčio pokytis tarp prognozuotų ir eksperimentinių\n struktūrų didėjant struktūrų skirtumui tarp eksperimentinių baltymų struktūrų")


print(figure5d_filtered_min)
ggsave("Stats/fig5d_full_min_a3_RMSD_RMSD.png", plot = figure5d_filtered_min, width = 8, height = 6)


# p-value = 0.4007, cor = 0.09845275 
cor.test(filtered_data$APO_vs_HOLO, filtered_data$Max_RMSD,
         method = 'pearson')

# p-value = 0.7313, cor = -0.04031271 
cor.test(filtered_data$APO_vs_HOLO, filtered_data$Min_RMSD,
         method = 'pearson')

combined_figure = grid.arrange(
  figure5a_full_max, figure5b_full_min, figure5c_filtered_max, figure5d_filtered_min,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig5a_combined.png", plot = combined_figure, width = 16, height = 8)


# FIG 5A

global_rmsd <- read.csv("reviewed_RMSD_info.csv", sep=',')  # Replace with your actual file path
global_rmsd$holo_id <- tolower(global_rmsd$holo_id)
data_with_global_rmsd <- merge(
  data,
  global_rmsd,
  by.x  = "Holo_PDB_ID",        # key in df1_sub
  by.y  = "holo_id",  # key in df2_sub
  all.x = TRUE          # left join (keeps all rows of df1_sub)
)
data_with_global_rmsd$max_cluster_RMSD = pmax(data_with_global_rmsd$holo_RMSD_max, data_with_global_rmsd$apo_RMSD_max)
data_with_global_rmsd$min_cluster_RMSD = pmin(data_with_global_rmsd$holo_RMSD_min, data_with_global_rmsd$apo_RMSD_min)

data_with_global_rmsd$min_cluster_RMSD = pmin(data_with_global_rmsd$holo_RMSD_min, data_with_global_rmsd$apo_RMSD_min)


# data_with_global_rmsd <- data_with_global_rmsd %>%
#   rowwise() %>%
#   mutate(row_max = max(c_across(c(holo_RMSD_max, apo_RMSD_max)), na.rm = TRUE)) %>%
#   ungroup()
print(mean(data_with_global_rmsd$max_cluster_RMSD - data_with_global_rmsd$min_cluster_RMSD))

data_with_global_rmsd$family = ifelse((data_with_global_rmsd$max_cluster_RMSD - data_with_global_rmsd$min_cluster_RMSD) < 1.273556, "Homogeneous", "Heterogeneous")
table(data_with_global_rmsd$family)


# Create the box plot
figure6a <- ggplot(data_with_global_rmsd, aes(x = family, y = Max_RMSD)) +
  geom_boxplot() +
  labs(
    title = 'Maksimalių Cα-RMSD įverčių pasiskirstymas baltymų šeimose',
    # title = "Box Plot of lowest RMSD Values by Family Group",
    x = "Šeima",
    y = "Maksimalus RMSD įvertis"
  ) +
  scale_x_discrete(labels = c("Heterogeniški", "Homogeniški")) +
  scale_y_continuous(limits = c(0,5)) +
  theme_minimal()
print(figure6a)
ggsave("Stats/RMSD_vs_families_max.png", plot = figure6a, width = 8, height = 6)

data_homogeneous <- subset(data_with_global_rmsd, family == "Homogeneous")
data_heterogeneous <- subset(data_with_global_rmsd, family == "Heterogeneous")

print(mean(data_homogeneous$Max_RMSD))
print(mean(data_heterogeneous$Max_RMSD))

wilcox_test <- wilcox.test(
  data_homogeneous$Max_RMSD,
  data_heterogeneous$Max_RMSD,
  alternative = "less"
)


# Print the Wilcoxon test results
print("Wilcoxon Rank-Sum Test Results:")
print(wilcox_test)


# Create the box plot
figure6b <- ggplot(data_with_global_rmsd, aes(x = family, y = Min_RMSD)) +
  geom_boxplot() +
  labs(
    title = 'Minimalių Cα-RMSD įverčių pasiskirstymas baltymų šeimose',
    # title = "Box Plot of lowest RMSD Values by Family Group(Apo Alpha)",
    x = "Šeima",
    y = "Minimalus RMSD įvertis"
  ) +
  scale_x_discrete(labels = c("Heterogeniški", "Homogeniški")) +
  scale_y_continuous(limits = c(0,5)) +
  theme_minimal()

print(figure6b)
ggsave("Stats/RMSD_vs_families_min.png", plot = figure6b, width = 8, height = 6)

combined_figure = grid.arrange(
  figure6a, figure6b,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig6a_combined.png", plot = combined_figure, width = 16, height = 8)



wilcox_test <- wilcox.test(
  data_homogeneous$Min_RMSD,
  data_heterogeneous$Min_RMSD,
  alternative = "less"
)

# Print the Wilcoxon test results
print("Wilcoxon Rank-Sum Test Results:")
print(wilcox_test)

# max_vid(RMSD)
print(mean(data_with_global_rmsd$max_cluster_RMSD))
data_with_global_rmsd$flexibility = ifelse((data_with_global_rmsd$max_cluster_RMSD) < 1.678111, "rigid", "flexible")
table(data_with_global_rmsd$flexibility)

figure7a <- ggplot(data_with_global_rmsd, aes(x = flexibility, y = Max_RMSD)) +
  geom_boxplot() +
  labs(
    title = 'Maksimalių Cα-RMSD įverčių pasiskirstymas baltymų šeimose',
    # title = "Box Plot of lowest RMSD Values by Family Group",
    x = "Lankstumas",
    y = "Maksimalus RMSD įvertis \n (tarp eksperimentinės struktūros ir prognozuotos struktūros)"
  ) +
  scale_x_discrete(labels = c("Lankstūs", "Nelankstūs")) +
  scale_y_continuous(limits = c(0,5)) +
  theme_minimal()
print(figure7a)
ggsave("Stats/RMSD_vs_flexibility_max.png", plot = figure7a, width = 8, height = 6)

data_flexible <- subset(data_with_global_rmsd, flexibility == "flexible")
data_rigid <- subset(data_with_global_rmsd, flexibility == "rigid")

wilcox_test <- wilcox.test(
  data_rigid$Max_RMSD,
  data_flexible$Max_RMSD,
  alternative = "less"
)


# Print the Wilcoxon test results
print("Wilcoxon Rank-Sum Test Results:")
print(wilcox_test)

# Create the box plot
figure7b <- ggplot(data_with_global_rmsd, aes(x = flexibility, y = Min_RMSD)) +
  geom_boxplot() +
  labs(
    title = 'Minimalių Cα-RMSD įverčių pasiskirstymas lanksčiuose ir nelanksčiuose baltymuose',
    # title = "Box Plot of lowest RMSD Values by Family Group(Apo Alpha)",
    x = "Lankstumas",
    y = "Minimalus RMSD įvertis \n (tarp eksperimentinės struktūros ir prognozuotos struktūros)"
  ) +
  scale_x_discrete(labels = c("Lankstūs", "Nelankstūs")) +
  scale_y_continuous(limits = c(0,5)) +
  theme_minimal()

print(figure7b)
ggsave("Stats/RMSD_vs_flexibility_min.png", plot = figure7b, width = 8, height = 6)

combined_figure = grid.arrange(
  figure7a, figure7b,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig7a_combined.png", plot = combined_figure, width = 16, height = 8)



wilcox_test <- wilcox.test(
  data_rigid$Min_RMSD,
  data_flexible$Min_RMSD,
  alternative = "less"
)

# Print the Wilcoxon test results
print("Wilcoxon Rank-Sum Test Results:")
print(wilcox_test)



# — user parameters —
infile       <- "combined_apo_holo_rmsf.csv"
out_apo_png  <- "Stats/apo_rmsf_by_plddt.png"
out_holo_png <- "Stats/holo_rmsf_by_plddt.png"

# read data
df <- read_csv(infile)

# helper for bins
df <- df %>%
  filter(!is.na(plddt_apo)) %>%
  mutate(
    plddt_apo_grp  = cut(plddt_apo,  breaks = c(-Inf,50,70,90,Inf),
                         labels = c("≤ 50","50–70","70–90","> 90")),
    plddt_holo_grp = cut(plddt_holo, breaks = c(-Inf,50,70,90,Inf),
                         labels = c("≤ 50","50–70","70–90","> 90"))
  )

summary(df$plddt_apo)
table(is.na(df$plddt_apo))

# common theme
my_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(df$plddt_apo_grp)
# 1) Apo plot
figure8a <- ggplot(df, aes(x = plddt_apo_grp, y = rmsf)) +
  geom_boxplot() +
  # geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
  labs(
    title = "Apo Cα-RMSF suskirstyta į pLDDT grupes",
    x     = "pLDDT grupė",
    y     = "Cα-RMSF (Å)"
  ) +
  scale_y_continuous(limits = c(0,2.5)) +
  my_theme
print(figure8a)
ggsave(out_apo_png, figure8a, width = 6, height = 4, dpi = 300)

# 2) Holo plot
figure8b <- ggplot(df, aes(x = plddt_holo_grp, y = rmsf)) +
  geom_boxplot() +
  # geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
  labs(
    title = "Holo Cα-RMSF suskirstyta į pLDDT grupes",
    x     = "pLDDT grupė",
    y     = "Cα-RMSF (Å)"
  ) +
  scale_y_continuous(limits = c(0,2.5)) +
  my_theme
print(figure7b)
ggsave(out_holo_png, figure8b, width = 6, height = 4, dpi = 300)

combined_figure = grid.arrange(
  figure8a, figure8b,
  ncol = 2,
  widths = c(1,1)
)
ggsave("Stats/fig8a_combined.png", plot = combined_figure, width = 16, height = 8)



cat("Saved:\n",
    " - ", out_apo_png, "\n",
    " - ", out_holo_png, "\n", sep = "")

cor.test(df$plddt_apo, df$rmsf, method = "pearson")
cor.test(df$plddt_holo, df$rmsf, method = "pearson")

wilcox_test <- wilcox.test(
  df$plddt_apo,
  df$plddt_holo,
  alternative = "less"
)

# Print the Wilcoxon test results
print("Wilcoxon Rank-Sum Test Results:")
print(wilcox_test)
